import numpy as np
from abc import ABC, abstractmethod

"""
Object containing routines and data structures required to dynamically change the template matrix of a bundle based off a pre-determined
strategy
"""
class TempStrategy(ABC):

    def __init__(self, model, order=None):
        self.model = model
        self.ptope_hash = {}
        self.ptope_counter = 0
        self.order = order
        self.counter = 0
        
    """
    Method called before the transformation operation and maximization over parallotopes are performed.
    """
    @abstractmethod
    def open_strat(self, bund):
        pass

    """
    Method called after the transformation operation is performed.
    """
    @abstractmethod
    def close_strat(self, bund):
        pass

    def hash_ptope(self, temp_idxs, name=None):
        key = self.__generate_unique_id() if name is None else name
        self.ptope_hash[key] = temp_idxs
        return key

    def pop_temp(self, key):
        return self.ptope_hash.pop(key)

    def rm_ptope_from_bund(self, bund, key):
        ptope_labels = self.ptope_hash[key]
        bund.remove_dirs(self, ptope_labels)
        bund.remove_temp(self, key)

    def add_ptope_to_bund(self, bund,  key):
        ptope_labels = self.ptope_hash[key]
        #bund.add_dirs(self, dir_row_mat, ptope_labels)
        bund.add_temp(self, ptope_labels, key)

    def __generate_unique_id(self):
        self.ptope_counter += 1
        return "Ptope" + str(self.ptope_counter)
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

            strat.order = self.strat_freq[type(strat)]
            self.strat_list.append(strat)

    def open_strat(self, bund):

        for strat in self.strat_list:
            strat.open_strat(bund)

    def close_strat(self, bund):

        for strat in self.strat_list:
            strat.close_strat(bund)
