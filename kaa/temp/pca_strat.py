import numpy as np
from sklearn.decomposition import PCA

from kaa.templates import TempStrategy
from kaa.bundle import Bundle
from kaa.experiutil import generate_traj
from kaa.timer import Timer

"""
Implementation of creating templates through PCA
"""
class PCAStrat(TempStrategy):

    def __init__(self, model, traj_steps=5, num_trajs=100, iter_steps=10):
        super().__init__(model)

        self.traj_steps = traj_steps
        self.num_trajs = num_trajs
        self.iter_steps = iter_steps
        self.curr_temp = None


    def open_strat(self, bund):
        if not self.counter % self.iter_steps:
            print(f"OPEN TEMP: {bund.T}")

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

            dir_idxs = bund.add_dir(comps)

            self.hash_dir('PCADir',dir_idxs)
            self.curr_temp = dir_idxs
            
            #print("After:  L: {} \n T: {}".format(bund.L, bund.T))
            #print("After: offu: {} \n  offl: {}".format(bund.offu, bund.offl))

        return bund
        
    def close_strat(self, bund):

        if not self.counter % self.iter_steps:

            if self.counter:
                bund.remove_temp(self.temp_hash['PCATemp'])
                self.pop_temp('PCATemp')

            temp_idx = bund.add_temp([self.curr_temp])
            self.hash_temp('PCATemp', temp_idx)

        print(f"CLOSE TEMP MAT: {bund.T}")

        self.counter += 1
        return bund


    def __str__(self):
        return "PCA(Steps: {})".format(self.iter_steps)
