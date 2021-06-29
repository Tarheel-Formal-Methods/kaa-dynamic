import os
import csv

from examples.test_ll import test_arch_LL
from examples.test_cvdp import test_arch_CVDP
from examples.test_lovo import test_arch_LOVO21
from settings import KaaSettings

def run_all_benchmarks():
    csv_path = os.path.join(KaaSettings.DataDir, 'results.csv')
    with open(csv_path, 'w') as csv_file:
        writer = csv.writer(csv_file)

        print("Running CVDP20")
        cvdp_data = test_arch_CVDP()
        writer.writerow(['Kaa', 'CVDP20', 'mu1', 0, cvdp_data[0][0], cvdp_data[0][1]])
        writer.writerow(['Kaa', 'CVDP20', 'mu2', 0, cvdp_data[1][0], cvdp_data[1][1]])

        print("Running LAL020")
        ll_data = test_arch_LL()
        writer.writerow(['Kaa', 'LALO20', 'W001', 0, ll_data[0][0], ll_data[0][1][3]])
        writer.writerow(['Kaa', 'LALO20', 'W005', 0, ll_data[1][0], ll_data[1][1][3]])
        writer.writerow(['Kaa', 'LAL020', 'W01', 0, ll_data[2][0], ll_data[2][1][3]])

        print("Running LOVO21")
        lovo_data = test_arch_LOVO21()
        writer.writerow(['Kaa', 'L0V021', '', 0, lovo_data[0][0], lovo_data[0][1]])

if __name__ == '__main__':
    run_all_benchmarks()