import os
import csv

from examples.test_ll import test_arch_LL
from examples.test_lovo import test_arch_LOVO21
from settings import KaaSettings

def run_all_benchmarks():
    csv_path = os.path.join(KaaSettings.DataDir, 'results.csv')
    with open(csv_path, 'w') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(['Tool','Benchmark','Verified','Runtime', 'Final Width(LALO20)/Box Area (LOVO21)'])

        print("Running LAL020")
        ll_data = test_arch_LL()

        w001_time = ll_data[0][0].formatted_tot_time
        w005_time = ll_data[1][0].formatted_tot_time
        w01_time = ll_data[2][0].formatted_tot_time

        writer.writerow(['Kaa', 'LALO20', 'W001', 0, f"{w001_time[0]} min, {w001_time[1]} sec" , ll_data[0][1][3]])
        writer.writerow(['Kaa', 'LALO20', 'W005', 0, f"{w005_time[0]} min, {w005_time[1]} sec", ll_data[1][1][3]])
        writer.writerow(['Kaa', 'LAL020', 'W01', 0,  f"{w01_time[0]} min, {w01_time[1]} sec", ll_data[2][1][3]])

        print("Running LOVO21")
        lovo_data = test_arch_LOVO21()

        lovo_time = lovo_data[0][0].formatted_tot_time
        writer.writerow(['Kaa', 'L0V021', '', 0, f"{lovo_time[0]} min, {lovo_time[1]} sec" , lovo_data[0][1]])

if __name__ == '__main__':
    run_all_benchmarks()