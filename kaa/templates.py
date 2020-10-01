import numpy as np
from abc import ABC, abstractmethod

"""
Object containing routines and data structures required to dynamically change the template matrix of a bundle based off a pre-determined
strategy
"""
class TempStrategy(ABC):

    def __init__(self, model):
        self.model = model
        self.dir_hash = {}
        self.temp_hash = {}
        self.counter = 0
        
    """
    Method called before the transformation operation is made.
    """
    @abstractmethod
    def open_strat(self, bund):
        pass

    """
    Method called after the transformation operation is made.
    """
    @abstractmethod
    def close_strat(self, bund):
        pass

    def hash_dir(self, key, dir_idxs):
        self.dir_hash[key] = dir_idxs

    def pop_dir(self, key):
        return self.dir_hash.pop(key)

    def hash_temp(self, key, temp_idxs):
        self.temp_hash[key] = temp_idxs

    def pop_temp(self, key):
        return self.temp_hash.pop(key)



"""
This would just be the static strategy where we do not touch any of the bundles after initializing them.
"""
class StaticStrat(TempStrategy):

    def __init__(self, model):
        super().__init__(model)

    def open_strat(self, bund):
        return bund

    def close_strat(self, bund):
        return bund
