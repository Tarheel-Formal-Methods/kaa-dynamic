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

    def __init__(self, model, traj_steps=5, num_traj=100, iter_steps=10):
        super().__init__(model)
        self.dim = len(self.model.vars)
        self.counter = 0

        self.traj_steps = traj_steps
        self.num_traj = num_traj
        self.iter_steps = iter_steps


    def open_strat(self, bund):
        if not self.counter % self.iter_steps:

            if self.counter:
               'We remove the last template added through PCA'
               bund.remove_temp(-1)
               bund.remove_dir(self.dir_hash['PCADir'])
               self.pop_hash_dir('PCADir')

            #print("Before: L: {} \n T: {}".format(bund.L, bund.T))
            #print("Before: offu: {}  offl: {}".format(bund.offu, bund.offl))

            trajs = generate_traj(bund, self.num_traj, self.traj_steps)
            traj_mat = np.empty((self.num_traj, self.dim))
            for traj_idx, traj in enumerate(trajs):
                traj_mat[traj_idx] = traj.end_point

            pca = PCA(n_components=self.dim)
            pca.fit(traj_mat)

            comps = pca.components_
            mean = pca.mean_

            dir_idxs = bund.add_dir(comps)
            bund.add_temp([dir_idxs])
            self.hash_dir('PCADir',dir_idxs)
            #print("After:  L: {} \n T: {}".format(bund.L, bund.T))
            #print("After: offu: {} \n  offl: {}".format(bund.offu, bund.offl))
            
        self.counter += 1
        return bund
        
    def close_strat(self, bund):
        return bund


    def __str__(self):
        return "PCA(Steps: {})".format(self.iter_steps)
