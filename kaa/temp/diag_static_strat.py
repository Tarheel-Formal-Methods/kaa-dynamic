import numpy as np
from itertools import product

from kaa.templates import TempStrategy
from kaa.log import Output


class RandomDiagStaticStrat(TempStrategy):

    def __init__(self, model, num_diag_temps):
        super().__init__(model)
        self.num_diag_temps = num_diag_temps

    def open_strat(self, bund, step_num):
        if not step_num:
            temp_dirs = self.__gen_diag_dirs()

            for temp_num, temp_dirs_mat in enumerate(temp_dirs):
                temp_vec_labels = [f"RanDiagDir{temp_num}Dir{i}" for i in range(self.dim)]
                self.add_ptope_to_bund(bund, temp_dirs_mat, temp_vec_labels)

    def close_strat(self, bund, step_num):
        pass

    def reset(self):
        pass

    def __gen_diag_dirs(self):
        axis_dirs = np.zeros((self.dim, self.dim))

        for i in range(self.dim):
            axis_dirs[i][i] = 1

        'Generate all the positive, negative diags'
        diag_mat = np.copy(axis_dirs)
        for row_idx, row in enumerate(axis_dirs):
            prod_mat_num_rows = 2*(self.dim - (row_idx+1))

            if prod_mat_num_rows:
                prod_mat = np.empty((prod_mat_num_rows, self.dim))

                for prod_row_idx, prod_rows in enumerate(product([row], axis_dirs[row_idx+1:])):
                    prod_mat[2*prod_row_idx] = prod_rows[0] + prod_rows[1]
                    prod_mat[2*prod_row_idx+1] = prod_rows[0] - prod_rows[1]

                diag_mat = np.append(diag_mat, prod_mat, axis=0)

        temp_dirs = np.empty((self.num_diag_temps, self.dim, self.dim))

        for trial_num in range(self.num_diag_temps):
            cond_num = np.inf
            used_dir_indices = np.empty((self.num_diag_temps, self.dim))
            while cond_num > 10 or (ran_idxs.tolist() in used_dir_indices.tolist()):
                ran_idxs = np.random.choice(len(diag_mat), self.dim, replace=False)
                temp_mat = diag_mat[ran_idxs]
                cond_num = np.linalg.cond(temp_mat)
            #Output.prominent(f"Cond num calculated for Ptope {trial_num}: {cond_num}")
            used_dir_indices[trial_num] = ran_idxs

            temp_dirs[trial_num] = temp_mat

        return temp_dirs



