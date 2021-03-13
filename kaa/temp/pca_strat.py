import numpy as np
from sklearn.decomposition import PCA

from kaa.templates import TempStrategy, GeneratedDirs, SampledTrajData
from kaa.bundle import Bundle
from kaa.timer import Timer
from kaa.log import Output

"""
Abstract PCA class containing all of the tools PCA strats need.
"""
class AbstractPCAStrat(TempStrategy):

    def __init__(self, model, num_trajs):
        super().__init__(model)
        self.num_trajs = num_trajs

    def open_strat(self, bund, step_num):
        pass

    def close_strat(self, bund, step_num):
        pass

    """
    Generate PCA directions from Bundle object or fetch them from
    GeneratedPCADirs object.
    @params bund: Bundle object
            step_num: step number of computation
    @returns tuple with first component being the matrix containing the PCA directions
             and the second component being the corresponding labels to identify each direction generated.
    """
    def generate_pca_dir(self, bund, step_num):

        if self.dirs is None:
            #Output.bold_write(f"Generating Directons for PCA")
            trajs = self.generate_trajs(bund, self.num_trajs)
            traj_mat = trajs.end_points
            Output.bold_write(f"Done generating directons for PCA")


            #Output.bold_write(f"Running PCA")
            pca = PCA(n_components=self.dim) #Run PCA
            pca.fit(traj_mat)
            pca_dirs_mat = pca.components_
            #Output.bold_write(f"Done with PCA")

        else:
            #Output.bold_write(f"Fetching Directions")
            pca_dirs_mat = self.dirs.get_dirs_at_step(step_num) #Else fetch pre-generated directions.
            #Output.bold_write(f"Done fetching Directions")

        ptope_dir_labels = [str((step_num, comp_idx)) for comp_idx, _ in enumerate(pca_dirs_mat)]
        return pca_dirs_mat, ptope_dir_labels

"""
Implementation of creating templates through PCA
"""
class PCAStrat(AbstractPCAStrat):

    def __init__(self, model, num_trajs=300, iter_steps=1):
        super().__init__(model, num_trajs)
        self.iter_steps = iter_steps
        self.ptope_queue = []

    """
    Opening PCA routine
    """
    def open_strat(self, bund, step_num):
        if not step_num % self.iter_steps:
            pca_comps, pca_comp_labels  = self.generate_pca_dir(bund, step_num)

            'Add the components to the bundle and save last generated ptope data.'
            ptope_label = self.add_ptope_to_bund(bund, pca_comps, pca_comp_labels)
            self.ptope_queue.append(ptope_label)

    """
    Closing PCA routine
    """
    def close_strat(self, bund, step_num):
        if not step_num % self.iter_steps and step_num > 0:
            self.rm_ptope_from_bund(bund, self.ptope_queue.pop(0))

    """
    Reset the strategy for a new round of computation.
    """
    def reset(self):
        self.ptope_queue = []

    def __str__(self):
        return f"PCAStrat(Steps:{self.iter_steps})" if self.strat_order is None else f"PCAStrat{self.strat_order}(Steps:{self.iter_steps})"

"""
Delayed PCA
"""
class SlidingPCAStrat(AbstractPCAStrat):

    def __init__(self, model, num_steps=70, num_trajs=300, lifespan=3):
        super().__init__(model, num_trajs)

        self.pca_ptope_life_counter = {}
        self.lifespan = lifespan

    def open_strat(self, bund, step_num):
        #Output.bold_write(f"Calling Open Strat for {str(self)}")
        self.__add_new_ptope(bund, step_num)

        'Remove dead templates'
        for ptope_label in list(self.pca_ptope_life_counter.keys()):
            self.pca_ptope_life_counter[ptope_label] -= 1

            if self.pca_ptope_life_counter[ptope_label] == 0:
                self.rm_ptope_from_bund(bund, ptope_label)
                self.pca_ptope_life_counter.pop(ptope_label)

        assert len(self.pca_ptope_life_counter) == min(step_num + 1, self.lifespan), f"Number of templates don't match. {len(self.pca_ptope_life_counter)} at step num {step_num}"

    def close_strat(self, bund, step_num):
        #Output.bold_write(f"Calling Closing Strat for {str(self)}")
        pass

    """
    Auxiliary method to generate PCA directions and add them as
    directions for a new ptope every step.
    @params bund: Bundle object
            step_num: step number in computation
    @returns: label string for generated ptope
    """
    def __add_new_ptope(self, bund, step_num):
        new_pca_dirs, new_dir_labels = self.generate_pca_dir(bund, step_num)
        new_ptope_label = self.add_ptope_to_bund(bund, new_pca_dirs, new_dir_labels)
        self.pca_ptope_life_counter[new_ptope_label] = self.lifespan + 1#Add fresh ptope and lifespan to step list

        return new_ptope_label

    """
    Reset the strategy for a new round of computation.
    """
    def reset(self):
         self.pca_ptope_life_counter = {}

    def __str__(self):
        return f"SlidingPCAStrat(Lifespan:{self.lifespan})"

class GeneratedPCADirs(GeneratedDirs):

    def __init__(self, model, num_steps, num_trajs, dir_mat=None, sampled_points=None):
        self.num_trajs = num_trajs
        self.num_steps = num_steps

        if dir_mat is None or sampled_points is None:
            pca_dirs, traj_points = self.__generate_pca_dir(model)
            super().__init__(model, pca_dirs, traj_points)
        else:
            super().__init__(model, dir_mat, sampled_points)

    """
    Pre-generates all PCA directions using random initial points from the initial box.
    @params model: Model object containing dynamics and initial box
    @returns matrix containing pre-generated PCA dirs.
    """
    def __generate_pca_dir(self, model):
        bund = model.bund
        dim = model.dim

        generated_pca_dir_mat = np.empty((dim*self.num_steps, dim))
        trajs = bund.getIntersect().generate_ran_trajs(self.num_trajs, self.num_steps) #trajs is TrajCollecton object'

        for step in range(self.num_steps):
            pca = PCA(n_components=dim)
            pca.fit(trajs[step]) #Takes point data from the (num_step)-th step of trajectories contained in TrajCollecton

            pca_dirs = pca.components_
            generated_pca_dir_mat[step*dim : (step + 1)*dim,] = pca_dirs

        return generated_pca_dir_mat, trajs.get_mat()
