import numpy as np
from itertools import product
from random import uniform
from sklearn.decomposition import PCA

from kaa.templates import TempStrategy
from kaa.bundle import Bundle
from kaa.experiutil import generate_traj
from kaa.timer import Timer

"""
WARNING: FRANKSTEIN'D CODE FROM PCA/LINAPP STRAT. ONLY MEANT TO BE A TEMPORARY FIX UNTIL MORE VERSATIBLE/EXTENIBLE SOLUTION IS DESIGNED.
"""
class PCALinStrat(TempStrategy):

    def __init__(self, model, traj_steps=5, num_trajs=100, iter_steps=10, cond_threshold=15):
        super().__init__(model)

        self.traj_steps = traj_steps
        self.num_trajs = num_trajs
        self.iter_steps = iter_steps
        self.dim = model.dim
        
        self.unit_dir_mat = np.zeros((self.dim, self.dim))
        for i in range(self.dim):
            self.unit_dir_mat[i][i] = 1

        self.cond_threshold = cond_threshold

    def open_strat(self, bund):
        if not self.counter % self.iter_steps:
            #print(f"OPEN TEMP: {bund.T}")

            #print("Before: L: {} \n T: {}".format(bund.L, bund.T))
            #
            #print("Before: offu: {}  offl: {}".format(bund.offu, bund.offl))

            trajs = generate_traj(bund, self.num_trajs, self.traj_steps)
            traj_mat = np.empty((self.num_trajs, bund.dim))

            'Populate data matrix for PCA routine.'
            for traj_idx, traj in enumerate(trajs):
                traj_mat[traj_idx] = traj.end_point

            pca = PCA(n_components=self.model.dim)
            pca.fit(traj_mat)

            comps = pca.components_
            mean = pca.mean_

            approx_A = self._approx_A(bund, self.dim+2)
            inv_A = np.linalg.inv(approx_A)
            lin_dir = np.dot(self.unit_dir_mat, inv_A)
            """
            cond_num = np.linalg.cond(lin_dir)
            print(f"COND NUM: {cond_num}")

            if cond_num > self.cond_threshold:
                closest_dirs = self.__find_closest_dirs(lin_dir)
                lin_dir = self.__merge_closest_dirs(lin_dir, closest_dirs)
                print(f"LIN DIR: {lin_dir}")
            """
            lin_dir_idxs = bund.add_dir(lin_dir)
            self.unit_dir_mat = lin_dir

            pca_dir_idxs = bund.add_dir(comps)

            if not self.counter:
                bund.add_temp([pca_dir_idxs, lin_dir_idxs])
                self.hash_dir('PrevDir', pca_dir_idxs +  lin_dir_idxs)

            #print("After:  L: {} \n T: {}".format(bund.L, bund.T))
            #print("After: offu: {} \n  offl: {}".format(bund.offu, bund.offl))

        return bund

    def close_strat(self, bund):

        if not self.counter % self.iter_steps:
            if self.counter:
                bund.remove_dir(self.dir_hash['PrevDir'])

        #print(f"CLOSE TEMP MAT: {bund.T}")

        self.counter += 1
        return bund



    def _approx_A(self, bund, num_traj):

        trajs = generate_traj(bund, num_traj, self.iter_steps)
        coeff_mat = np.zeros((self.dim*num_traj,self.dim**2), dtype='float')

        'Initialize the A matrix containing linear constraints.'
        for t_idx, t in enumerate(trajs):
            for i in range(self.dim):
                for j in range(self.dim):
                    coeff_mat[i+self.dim*t_idx][i*self.dim+j] = t.start_point[j]

        b_mat_vec = np.asarray([t.end_point for t in trajs], dtype='float')
        b_mat = b_mat_vec.flatten()

        m = np.linalg.lstsq(coeff_mat, b_mat, rcond=None)[0]
        return m.reshape((self.dim,self.dim))

    def __merge_closest_dirs(self, dir_mat, closest_dirs):
        first_dir, second_dir = (0,1)
        merged_dir = (dir_mat[first_dir] + dir_mat[second_dir]) / 2
        ortho_dir = [ uniform(-1,1) for _ in range(self.dim) ]

        accum_mat = np.delete(dir_mat, [0,1], axis=0)
        accum_mat = np.append(accum_mat, [merged_dir], axis=0)


        'Orthogonalize generated vector through Gram-Schmidt'
        for row in accum_mat:
            ortho_dir -= np.dot(ortho_dir, row) * row

        norm_ortho_dir = ortho_dir / np.linalg.norm(ortho_dir)
        return np.append(accum_mat, [norm_ortho_dir], axis=0)

    def __find_closest_dirs(self, dir_mat):

        closest_pair = None
        closest_dot_prod = 0

        norm_dir_mat = self.__normalize_mat(dir_mat)

        print(f"NORM_DIR_MAT: {norm_dir_mat}")

        for (first_idx, first_dir), (second_idx, second_dir) in product(zip(range(len(norm_dir_mat)), norm_dir_mat), repeat=2):
            curr_product = np.dot(first_dir, second_dir)

            #print(curr_product)
            if curr_product > closest_dot_prod and first_idx != second_idx:
                closest_pair = (first_idx, second_idx)
                closest_dot_prod = curr_product

        return closest_pair

    def __normalize_mat(self, mat):
        return mat / np.linalg.norm(mat, ord=2, axis=1, keepdims=True)



    def __str__(self):
        return "LinApp(Steps:{})".format(self.iter_steps)
