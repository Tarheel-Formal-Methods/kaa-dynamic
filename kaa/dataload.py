import os
from datetime import datetime
from collections import OrderedDict
import csv

from kaa.settings import KaaSettings


class CovidDataLoader:

    def __init__(self, filename):
        self.csv_path = os.path.join(KaaSettings.DataDir, filename)

    def __load_data(self, earlier_date, later_date):
        data_dict = OrderedDict()

        with open(self.csv_path) as file:
            reader = csv.reader(file, delimiter=',')
            for line, row in enumerate(reader):
                if line > 0:
                    row_time_obj = datetime.strptime(row[0], "%d %B %Y")

                    if earlier_date <= row_time_obj <= later_date:
                        converted_date = row_time_obj.strftime("%D")
                        data_dict[converted_date] = OrderedDict()

                        data_dict[converted_date]['confirmed'] = row[3]
                        data_dict[converted_date]['recovered'] = row[5]
                        data_dict[converted_date]['deceased'] = row[7]

        return data_dict

    def fetch_data(self, earlier_date, later_date):
        earlier_date_obj = datetime.strptime(earlier_date, "%m/%d/%y") if earlier_date else None
        later_date_obj = datetime.strptime(later_date, "%m/%d/%y") if later_date else None

        if earlier_date and later_date:
            assert earlier_date_obj < later_date_obj, "Later date must be strictly in future"

        return self.__load_data(earlier_date_obj, later_date_obj)
