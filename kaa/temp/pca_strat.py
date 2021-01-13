import numpy as np
from sklearn.decomposition import PCA

from kaa.templates import TempStrategy, GeneratedDirs
from kaa.bundle import Bundle
from kaa.timer import Timer

"""
Abstract PCA class containing all of the tools PCA strats need.
"""
class AbstractPCAStrat(TempStrategy):

    def __init__(self, model, num_trajs, pca_dirs):
        assert pca_dirs is None or isinstance(pca_dirs, GeneratedPCADirs), "PCA Strategies may only take pre-generated PCA directions."
        super().__init__(model)
        self.num_trajs = num_trajs
        self.pca_dirs = pca_dirs

    def open_strat(self, bund, step_num):
        pass

    def close_strat(self, bund, step_num):
        pass

    def generate_pca_dir(self, bund, step_num):
        if self.pca_dirs is None:
            trajs = bund.getIntersect().generate_traj(self.num_trajs, 1, sample=False)
            traj_mat = trajs.end_points

            pca = PCA(n_components=self.dim)
            pca.fit(traj_mat)
            pca_dirs_mat = pca.components_

        else:
            pca_dirs_mat = self.pca_dirs.get_dirs_at_step(step_num + 1)

        ptope_dir_labels = [str((step_num, comp_idx)) for comp_idx, _ in enumerate(pca_dirs_mat)]
        return pca_dirs_mat, ptope_dir_labels

"""
Implementation of creating templates through PCA
"""
class PCAStrat(AbstractPCAStrat):

    def __init__(self, model, num_trajs=300, iter_steps=1, pca_dirs=None):
        super().__init__(model, num_trajs, pca_dirs)
        self.iter_steps = iter_steps
        self.last_ptope = None

    def open_strat(self, bund, step_num):
        if not step_num % self.iter_steps:
            #print(f"OPEN DIR/TEMP MAT: {bund.L},  {bund.T}")
            pca_comps, pca_comp_labels  = self.generate_pca_dir(bund, step_num)

            'Add the components to the bundle.'
            ptope_label = self.add_ptope_to_bund(bund, pca_comps, pca_comp_labels)
            self.last_ptope = ptope_label

    def close_strat(self, bund, step_num):
        if not step_num % self.iter_steps and step_num > 0:
                self.rm_ptope_from_bund(bund, self.last_ptope)

    def reset(self):
        self.last_ptope = None

    def __str__(self):
        return f"PCAStrat(Steps:{self.iter_steps})" if self.strat_order is None else f"PCAStrat{self.strat_order}(Steps:{self.iter_steps})"

"""
Delayed PCA
"""
class SlidingPCAStrat(AbstractPCAStrat):

    def __init__(self, model, num_steps=70, num_trajs=300, lifespan=3, pca_dirs=None):
        super().__init__(model, num_trajs, pca_dirs)
        self.pca_ptope_life_data = {}
        self.life_span = lifespan

    def open_strat(self, bund, step_num):
        self.__add_new_ptope(bund, step_num)
        #print("Before:  L: {} \n T: {}".format(bund.L, bund.T))
        #print("Before: offu: {} \n  offl: {}".format(bund.offu, bund.offl))

    def close_strat(self, bund, step_num):
        'Remove dead templates'
        for ptope_label in list(self.pca_ptope_life_data.keys()):
            self.pca_ptope_life_data[ptope_label] -= 1
            if self.pca_ptope_life_data[ptope_label] == 0:
                self.rm_ptope_from_bund(bund, ptope_label)
                self.pca_ptope_life_data.pop(ptope_label)

        #print("After:  L: {} \n T: {}".format(bund.L, bund.T))
        #print("After: offu: {} \n  offl: {}".format(bund.offu, bund.offl))

    def __add_new_ptope(self, bund, step_num):
        new_pca_dirs, new_dir_labels = self.generate_pca_dir(bund, step_num)
        new_ptope_label = self.add_ptope_to_bund(bund, new_pca_dirs, new_dir_labels)
        self.pca_ptope_life_data[new_ptope_label] = self.life_span #Add fresh ptope and lifespan to step list

        return new_ptope_label

    def reset(self):
         self.pca_ptope_life_data = {}

    def __str__(self):
        return f"SlidingPCAStrat(Lifespan:{self.life_span})" if self.strat_order is None else f"SlidingPCAStrat{self.strat_order}(Lifespan:{self.life_span})"

class GeneratedPCADirs(GeneratedDirs):

    def __init__(self, model, num_steps, num_trajs):
        super().__init__(model, self.__generate_pca_dir(model, num_trajs, num_steps))

    def __generate_pca_dir(self, model, num_trajs, num_steps):
        bund = model.bund
        dim = model.dim

        generated_pca_dir_mat = np.empty((dim*num_steps, dim))
        trajs = bund.getIntersect().generate_traj(num_trajs, num_steps, sample=False) #trajs is TrajCollecton object'

        for step in range(num_steps):
            pca = PCA(n_components=dim)
            pca.fit(trajs[step]) #Takes point data from the (num_step)-th step of trajectories contained in TrajCollecton

            pca_dirs = pca.components_
            generated_pca_dir_mat[step*dim:(step+1)*dim,] = pca_dirs

        return generated_pca_dir_mat
