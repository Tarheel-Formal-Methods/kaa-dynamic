import numpy as np
from abc import ABC, abstractmethod
from collections import namedtuple
from kaa.settings import KaaSettings

class SampledTrajData:

    def __init__(self, initial_points, image_points):
        self.initial_points = initial_points
        self.image_points = image_points

    def __getitem__(self, idx):
        TrajData = nampledtuple('TrajData', ['initial_points, image_points'])
        return TrajData(self.initial_points[idx], self.image_points[idx])

"""
Object containing routines and data structures required to dynamically change the template matrix of a bundle based off a pre-determined
strategy.
The ptope dictionary represents the local knowledge the strategy
has on the ptopes which it creates and modifies.
"""
class TempStrategy(ABC):

    def __init__(self, model, use_supp_points=False, stratorder=None):
        self.model = model
        self.dim = model.dim
        self.ptope_hash = {}
        self.ptope_counter = 0
        self.strat_order = stratorder
        self.initial_points = []
        self.image_points = []
        self.use_supp_points = use_supp_points
        self.pregen_dirs = None
        self.dirs = None


    """
    Method called before the transformation operation and maximization over parallotopes are performed.
    """
    @abstractmethod
    def open_strat(self, bund, step_num):
        pass

    """
    Method called after the transformation operation is performed.
    """
    @abstractmethod
    def close_strat(self, bund, step_num):
        pass

    """
    Method to reset all local tabular attributes.
    """
    @abstractmethod
    def reset(self):
        pass

    """
    Inserts list of direction labels associated to a ptope into the ptope dictonary.
    The method returns a label for the ptope if the name is not specified.
    @params dir_label_list: list of labels to direction composing the template.
    @returns key of hashed ptope
    """
    def __hash_ptope(self, dir_label_list, name=None):
        assert len(dir_label_list) == self.dim, "Number of directions associated to a parallelotope must match the dimension of the system."
        
        key = self.__generate_unique_key() if name is None else name
        self.ptope_hash[key] = dir_label_list
        return key

    """
    Removes an entry in the ptope dictionary with specified key.
    @params key: key of ptope to be removed
    """
    def __pop_ptope(self, key):
        return self.ptope_hash.pop(key)

    """
    Removes a template contained in the ptope dictionary from an input bundle
    @params bund: input Bundle object
            key: key of ptope indexed in self.ptope_hash
    """
    def rm_ptope_from_bund(self, bund, ptope_label):
        ptope_dir_labels = self.ptope_hash[ptope_label]
        bund.remove_dirs(self, ptope_dir_labels)
        bund.remove_temp(self, ptope_label)
        self.__pop_ptope(ptope_label)

    """
    Adds a template contained in the ptope dictionary from an input bundle
    @params bund: input Bundle object
            key: key of ptope indexed in self.ptope_hash
    """
    def add_ptope_to_bund(self, bund,  dir_row_mat, dir_labels):
        ptope_name = self.__hash_ptope(dir_labels)
        bund.add_dirs(self, dir_row_mat, dir_labels)
        bund.add_temp(self, dir_labels, ptope_name)

        return ptope_name
        
    """
    Generates a unique key for a ptope.
    @returns unique label for ptope
    """
    def __generate_unique_key(self):
        self.ptope_counter += 1
        return "Ptope" + str(self.ptope_counter)
    
    """
    Generate trajectories according to scheme defined in KaaSettings
    @params bund: Bundle object
            num_trajs: number of trajectories to generate
    @returns TrajCollection object wrapping generated Traj objects.
    """
    def generate_trajs(self, bund, num_trajs):
        
        if self.use_supp_points:
            trajs = bund.getIntersect().generate_supp_trajs(bund.L, 1)
        else:
            trajs = bund.getIntersect().generate_ran_trajs(num_trajs, 1)

        self.initial_points.append(trajs.start_points)
        self.image_points.append(trajs.end_points)
        return trajs
    
    """
    Fetch data containing the initial points of the trajectories used and the points where they were propagated.
    @returns TrajData object containing all relevant data
    """
    def fetch_traj_data(self):
        
        if self.dirs is None:
            return SampledTrajData(self.initial_points, self.image_points)
        else:
            sampled_traj_pts = self.dirs.sampled_points
            #print(f"Shape of Sampled Points: {self.dirs.sampled_points.shape}")
            return SampledTrajData(sampled_traj_pts,
                                   sampled_traj_pts[:,1:,:]) #Edge cases here to consider.

"""
This would just be the static strategy where we do not touch any of the bundles after initializing them.
"""
class StaticStrat(TempStrategy):

    def __init__(self, model):
        super().__init__(model)

    def open_strat(self, bund, step_num):
        return bund

    def close_strat(self, bund, step_num):
        return bund

"""
A wrapper enveloping multiple strategies working in tandem.
"""
class MultiStrategy(TempStrategy):

    def __init__(self, *var_tup):
        self.strat_list = []
        self.strat_freq = {}

        for strat in var_tup:

            assert isinstance(strat, TempStrategy), "Only a list of TempStrategy objects must be passed into a MultiStrategy."

            if type(strat) in self.strat_freq:
                self.strat_freq[type(strat)] += 1
            else:
                self.strat_freq[type(strat)] = 1

            strat.strat_order = self.strat_freq[type(strat)]
            self.strat_list.append(strat)

    @property
    def strats(self):
        return self.strat_list

    def open_strat(self, bund, step_num):
        for strat in self.strat_list:
            strat.open_strat(bund, step_num)

    def close_strat(self, bund, step_num):
        for strat in self.strat_list:
            strat.close_strat(bund, step_num)

    def reset(self):
        for strat in self.strat_list:
            strat.reset()

    def __str__(self):
        return ' and '.join([str(strat) for strat in self.strat_list])

"""
Wrapper over matrix of pre-generated dirs.
"""
class GeneratedDirs:

    def __init__(self, model, dir_mat, sampled_points):
        self.model = model
        self.dim = model.dim
        self.dir_mat = dir_mat
        self.num_steps = len(dir_mat)
        self.sampled_points = sampled_points

    @classmethod
    def from_mat(cls, model, dir_mat, sampled_points):
        return cls(model, dir_mat, sampled_points)

    def get_dirs_at_step(self, step_num):
        return self.dir_mat[step_num*self.dim : (step_num + 1)*self.dim,]
