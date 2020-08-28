import numpy as np
from sklearn.decomposition import PCA

from kaa.templates import TempStrategy, Template, Direction
from kaa.bundle import Bundle
from kaa.experiutil import generate_traj

iter_steps = 10
num_traj = 30
traj_steps = 10

"""
Implementation of creating templates through PCA
"""
class PCAStrat(TempStrategy):

    def __init__(self, model):
        super().__init__(model)
        self.dim = len(self.model.vars)
        self.counter = 0

    def open_strat(self, bund):
        T = Template(bund)
        L = Direction(bund)
        if not self.counter % iter_steps:
            #TrajSet wrapper?
            #Better memory management of Template class
            #Hacky revise Matrix class plz and revise this class def.

            if self.counter:
               T.remove(-1)
               for _ in range(self.dim):
                   L.remove(-1)
                   bund.offu = np.delete(bund.offu, -1, axis=0)
                   bund.offl = np.delete(bund.offl, -1, axis=0)
                   
            print("L: {}  T: {}".format(L.M, T.M))
            print("offu: {}  offl: {}".format(bund.offu, bund.offl))

            trajs = generate_traj(bund, num_traj, traj_steps)
            traj_mat = np.empty((num_traj, self.dim))
            for traj_idx, traj in enumerate(trajs):
                traj_mat[traj_idx] = traj[-1]

            pca = PCA(n_components=self.dim)
            pca.fit(traj_mat)

            comps = pca.components_
            mean = pca.mean_

            #print(comps)
            dir_idx = [ L.add([c]) for c in comps ]
            T.add([dir_idx])

            print("L: {}  T: {}".format(L.M, T.M))

            offu = np.append(bund.offu, [0 for _ in range(self.dim)])
            offl = np.append(bund.offl, [0 for _ in range(self.dim)])

            print("offu: {}  offl: {}".format(offu, offl))

            b = Bundle(bund.model, T.M, L.M, offu, offl)

        else:
            b =  bund
            
        self.counter += 1
        return b
        
    def close_strat(self, bund):
        return bund
