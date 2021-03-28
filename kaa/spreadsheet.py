from openpyxl import Workbook
from openpyxl.worksheet.table import Table, TableStyleInfo
from pathlib import Path
from datetime import date
import os

from kaa.settings import PlotSettings

class SpreadSheetTable:

    def __init__(self, sheet, data_label, inputs, num_trials, offset):
        self.row_offset = offset
        self.column_offset = 0
        self.sheet = sheet
        self.inputs = inputs
        self.label = data_label
        self.num_trials = num_trials
        self.row_dict = None
        self.num_rows = 0
        self.num_cols = 0

    def generate_table(self):
        self.__generate_headers()

        """
        Initialize label-row dictionary in relative ordering. (0,0) top left most corner starts with entry
        corresponding to first experiment, first trial.
        """
        self.row_dict = {experi_input['label'] : row for row, experi_input in enumerate(self.inputs)}

        for experi_label, row in self.row_dict.items():
            header_coord = self.__xlsx_tab_coord(row + 1, 0)
            self.sheet[header_coord] = experi_label

            if self.num_trials > 1:
                self[row, self.num_trials] = f"=AVERAGE({self.__row_xlsx_statement(row)})"
                self[row, self.num_trials+1] = f"=STDEV({self.__row_xlsx_statement(row)})"

        data_table = Table(displayName=f"{self.label}Table", ref=self.__get_table_corners())
        style = TableStyleInfo(name="TableStyleMedium9",  showRowStripes=True)

        return data_table

    def save_data_into_table(self, experi_label, trial_num, data):
        row = self.row_dict[experi_label]
        self[row, trial_num] = data

    def __row_xlsx_statement(self, row):
        start_cell = self.__rel2abs(row, 0)
        end_cell = self.__rel2abs(row, self.num_trials-1)
        return f"{self.__xlsx_tab_coord(*start_cell)}:{self.__xlsx_tab_coord(*end_cell)}"

    def __generate_headers(self):
        headers = ["Strategy"] + [f"Trial {i+1}" for i in range(self.num_trials)]
        if self.num_trials > 1:
            headers += ["Mean", "Stdev"]

        for column, header in enumerate(headers):
            self.sheet[self.__xlsx_tab_coord(0, column)] = header

        self.num_rows = len(self.inputs) + 1
        self.num_cols = len(headers)

    def __rel2abs(self, x, y):
        return (x+1, y+1)

    def __xlsx_tab_coord(self, x, y):
        return  chr(66 + self.column_offset + y) + str(1 + self.row_offset + x)

    def __get_table_corners(self):
        return f"{self.__xlsx_tab_coord(0,0)}:{self.__xlsx_tab_coord(self.num_rows-1, self.num_cols-1)}"

    def __getitem__(self, key):
        abs_coord = self.__rel2abs(*key)
        return self.sheet[self.__xlsx_tab_coord(*abs_coord)]

    def __setitem__(self, key, value):
        abs_coord = self.__rel2abs(*key)
        self.sheet[self.__xlsx_tab_coord(*abs_coord)] = value

"""
Data class for storing objects relevant to saving data to openpyxl spreadsheets
"""
class SpreadSheet:

    def __init__(self, model, name):
        self.model = model
        self.name = name
        self.workbook = Workbook()
        self.sheet = self.workbook.active
        self.num_trials = None
        self.table_dict = {}
        self.data_pwd = self.__gen_data_directory()

    """
    Saves data into a desired cell in spreadsheet.
    """
    def save_data_into_sheet(self, experi_label, trial_num, data_tup):
        assert len(data_tup) == len(self.table_dict), "Data tuple must contain a data entry for each table in the spreadsheet."

        for label, data in zip(self.table_dict.keys(), data_tup):
            table = self.table_dict[label]
            table.save_data_into_table(experi_label, trial_num, data)

        self.workbook.save(filename=os.path.join(self.data_pwd, self.name + '.xlsx'))

    """
    Generates directory path used to save spreadsheet into disk.
    @returns total path to data directory
    """
    def __gen_data_directory(self):
        data_pwd = os.path.join(PlotSettings.default_fig_path, 'Spreadsheets', date.today().isoformat(), str(self.model))
        Path(data_pwd).mkdir(parents=True, exist_ok=True)
        return data_pwd

    """
    Initializes openpyxl spreadsheet to dump resulting data.
    """
    def generate_sheet(self, inputs, data_labels, num_trials):
        self.num_trials = num_trials

        for label_order, label in enumerate(data_labels):
            tab = SpreadSheetTable(self.sheet, label, inputs, num_trials, label_order*(len(inputs)+2))
            self.table_dict[label] = tab

            pxyl_table = tab.generate_table()
            self.sheet.add_table(pxyl_table)

        self.workbook.save(filename=os.path.join(self.data_pwd, self.name + '.xlsx'))


    def __xlsx_coord(self, x, y):
        return chr(66 + y) + str(x+1)

    def __getitem__(self, key):
        return self.sheet[self.__xlsx_coord(*key)]

    def __setitem__(self, key, value):
        self.sheet[self.__xlsx_coord(*key)] = value
