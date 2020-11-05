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
        self.pca_ptope_queue = []

    def open_strat(self, bund):

        if not self.counter % self.iter_steps:
            #print(f"OPEN DIR/TEMP MAT: {bund.L},  {bund.T}")

            trajs = generate_traj(bund, self.num_trajs, self.traj_steps)
            traj_mat = np.empty((self.num_trajs, bund.dim))

            'Populate data matrix for PCA routine.'
            for traj_idx, traj in enumerate(trajs):
                traj_mat[traj_idx] = traj.end_point

            pca = PCA(n_components = self.model.dim)
            pca.fit(traj_mat)

            comps = pca.components_
            mean = pca.mean_

            ptope_labels = [str((counter, comp_idx)) for comp_idx, _ in enumerate(comps) ]
            bund.add_dir(self, comps, ptope_labels)

            'Hash labels and insert into Queue'
            self.pca_ptope_queue.append(self.hash_ptope(ptope_labels))
            
            #print("After:  L: {} \n T: {}".format(bund.L, bund.T))
            #print("After: offu: {} \n  offl: {}".format(bund.offu, bund.offl))

        return bund
        
    def close_strat(self, bund):

        if not self.counter % self.iter_steps:
            if self.counter:
                self.rm_ptope_from_bund(bund, self.pca_ptope_queue.pop(0))

            self.add_ptope_to_bund(bund, self.pca_ptope_queue[0])

        self.counter += 1
        return bund


    def __str__(self):
        return "PCA(Steps: {})".format(self.iter_steps)
