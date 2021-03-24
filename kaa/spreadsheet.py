from openpyxl import Workbook
from pathlib import Path
import os


"""
Data class for storing objects relevant to saving data to openpyxl spreadsheets
"""
class SpreadSheet:

    def __init__(self, model, name):
        self.model = model
        self.name = name
        self.workbook = Workbook()
        self.row_dict = None
        self.data_pwd = self.__gen_data_directory()

    """
    Saves data into a desired cell in spreadsheet.
    """
    def save_data_into_sheet(self, trial_num, num_trials, flow_label, data):
        column_offset = trial_num
        row_offset = row_dict[flow_label]

        sheet = self.workbook.active
        sheet[chr(66 + column_offset) + str(row_offset)] = data

        if column_offset == num_trials - 1:
            sheet[chr(66 + num_trials) + str(row_offset)] = f"=AVERAGE(B{row_offset}:{chr(66 + num_trials - 1)}{row_offset})"
            sheet[chr(66 + num_trials + 1) + str(row_offset)] = f"=STDEV(B{row_offset}:{chr(66 + num_trials - 1)}{row_offset})"

        self.workbook.save(filename=os.path.join(self.data_pwd, self.name + '.xlsx'))

    """
    Generates directory path used to save spreadsheet into disk.
    @returns total path to data directory
    """
    def __gen_data_directory(self):
        data_pwd = os.path.join(PlotSettings.default_fig_path, date.today().isoformat(), str(self.model))
        Path(data_pwd).mkdir(parents=True, exist_ok=True)
        return data_pwd

    """
    Initializes openpyxl spreadsheet to dump resulting data.
    """
    def generate_sheet(self, inputs, num_trials):
        sheet = self.workbook.active
        sheet.append(["Strategy"] + [f"Trial {i+1}" for i in range(num_trials)] + ["Mean", "Stdev"])

        'Initialize label-row dictionary'
        self.row_dict = {experi_input['label'] : row_idx + 2 for row_idx, experi_input in enumerate(inputs)}

        for experi_input in inputs:
            flow_label = experi_input['label']
            row = self.row_dict[flow_label]
            sheet['A' + str(row)] = flow_label

        self.workbook.save(filename=os.path.join(self.data_pwd, self.name + '.xlsx'))
