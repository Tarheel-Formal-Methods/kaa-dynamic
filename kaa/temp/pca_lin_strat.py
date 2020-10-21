import numpy as np
from sklearn.decomposition import PCA

from kaa.templates import TempStrategy
from kaa.bundle import Bundle
from kaa.experiutil import generate_traj
from kaa.timer import Timer

"""
WARNING: FRANKSTEIN'D CODE FROM PCA/LINAPP STRAT. ONLY MEANT TO BE A TEMPORARY FIX UNTIL MORE VERSATIBLE/EXTENIBLE SOLUTION IS DESIGNED.
"""
class PCALinStrat(TempStrategy):

    def __init__(self, model, traj_steps=5, num_trajs=100, iter_steps=10):
        super().__init__(model)

        self.traj_steps = traj_steps
        self.num_trajs = num_trajs
        self.iter_steps = iter_steps
        self.dim = model.dim
        
        self.unit_dir_mat = np.zeros((self.dim, self.dim))
        for i in range(self.dim):
            self.unit_dir_mat[i][i] = 1


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

    def __str__(self):
        return "LinApp(Steps:{})".format(self.iter_steps)
