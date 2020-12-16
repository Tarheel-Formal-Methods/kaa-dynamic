import numpy as np
from sklearn.decomposition import PCA

from kaa.templates import TempStrategy, GeneratedDirs
from kaa.bundle import Bundle
from kaa.timer import Timer

"""
Abstract PCA class containing all of the tools PCA strats need.
"""
class AbstractPCAStrat(TempStrategy):

    def __init__(self, model, traj_steps, num_trajs, pca_dirs):
        assert isinstance(pca_dirs, GeneratedPCADirs), "PCA Strategies may only take pre-generated PCA directions."

        super().__init__(model)
        self.traj_steps = traj_steps
        self.num_trajs = num_trajs
        self.pca_dirs = pca_dirs

    def open_strat(self, bund, step_num):
        pass

    def close_strat(self, bund, step_num):
        pass

    def generate_pca_dir(self, bund, step_num):
        if self.pca_dirs is None:
            trajs = bund.getIntersect().generate_traj(self.num_trajs, self.traj_steps)
            traj_mat = trajs.end_points

            pca = PCA(n_components=self.dim)
            pca.fit(traj_mat)
            pca_dirs_mat = pca.components_

        else:
            pca_dirs_mat = self.pca_dirs.get_dirs_at_step(step_num+1)

        ptope_dir_labels = [str((step_num, comp_idx)) for comp_idx, _ in enumerate(pca_dirs_mat)]
        return pca_dirs_mat, ptope_dir_labels

"""
Implementation of creating templates through PCA
"""
class PCAStrat(AbstractPCAStrat):

    def __init__(self, model, traj_steps=5, num_trajs=100, iter_steps=1, pca_dirs=None):
        super().__init__(model, traj_steps, num_trajs, pca_dirs)
        self.iter_steps = iter_steps
        self.pca_ptope_queue = []

    def open_strat(self, bund, step_num):

        if not step_num % self.iter_steps:
            #print(f"OPEN DIR/TEMP MAT: {bund.L},  {bund.T}")
            pca_comps, pca_comp_labels  = self.generate_pca_dir(bund, step_num)

            'Add the components to the bundle.'
            ptope_label = self.add_ptope_to_bund(bund, pca_comps, pca_comp_labels)
            self.pca_ptope_queue.append(ptope_label)

    def close_strat(self, bund, step_num):
        if not step_num % self.iter_steps:
            if step_num:
                last_ptope = self.pca_ptope_queue.pop(0)
                self.rm_ptope_from_bund(bund, last_ptope)
            
    def __str__(self):
        return f"PCAStrat(Steps:{self.iter_steps})" if self.strat_order is None else f"PCAStrat{self.strat_order}(Steps:{self.iter_steps})"

"""
Delayed PCA
"""
class SlidingPCAStrat(AbstractPCAStrat):

    def __init__(self, model, traj_steps=5, num_trajs=100, lifespan=3, pca_dirs=None):
        super().__init__(model, traj_steps, num_trajs, pca_dirs)
        self.pca_ptope_life = []
        self.life_span = lifespan

    def open_strat(self, bund, step_num):
        self.__add_new_ptope(bund, step_num)
        #print("Before:  L: {} \n T: {}".format(bund.L, bund.T))
        #print("Before: offu: {} \n  offl: {}".format(bund.offu, bund.offl))

    def close_strat(self, bund, step_num):
        'Remove dead templates'
        for ptope_label, life in self.pca_ptope_life:
            if life == 0:
                self.rm_ptope_from_bund(bund, ptope_label)

        self.pca_ptope_life = [(label, life-1) for label, life in self.pca_ptope_life if life > 0]

        #print("After:  L: {} \n T: {}".format(bund.L, bund.T))
        #print("After: offu: {} \n  offl: {}".format(bund.offu, bund.offl))

    def __add_new_ptope(self, bund, step_num):
        new_pca_dirs, new_dir_labels = self.generate_pca_dir(bund, step_num)
        new_ptope_label = self.add_ptope_to_bund(bund, new_pca_dirs, new_dir_labels)
        self.pca_ptope_life.append((new_ptope_label, self.life_span)) #Add fresh ptope and lifespan to step list

        return new_ptope_label
        
    def __str__(self):
        return f"DelayedPCAStrat(Lifespan:{self.life_span})" if self.strat_order is None else f"DelayedPCAStrat{self.strat_order}(Lifespan:{self.life_span})"

class GeneratedPCADirs(GeneratedDirs):

    def __init__(self, model, num_trajs, num_steps):
        super().__init__(model, self.__generate_pca_dir(model, num_trajs, num_steps))

    def __generate_pca_dir(self, model, num_trajs, num_steps):
        bund = model.bund
        dim = model.dim

        generated_pca_dir_mat = np.empty((dim*num_steps, dim))
        trajs = bund.getIntersect().generate_traj(num_trajs, num_steps) #trajs is TrajCollecton object'

        for step in range(num_steps):
            pca = PCA(n_components=dim)
            pca.fit(trajs[num_steps]) #Takes point data from the (num_step)-th step of trajectories contained in TrajCollecton

            pca_dirs = pca.components_
            generated_pca_dir_mat[step*dim:(step+1)*dim,] = pca_dirs

        return generated_pca_dir_mat
