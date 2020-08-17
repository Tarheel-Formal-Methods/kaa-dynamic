import numpy as np
from abc import ABC, abstractmethod

"""
Object wrapping over template matrix to ease dynamic template loading and unloading.
"""
class Template:

    def __init__(self, T):

        self.T = T

    def add_temp(self, temp_row):

        self.T = np.append(self.T, temp_row, axis=0)
        return len(self.T) - 1

    def remove_temp(self, index):

        self.T = np.delete(self.T, index, axis=0)

"""
Object containing routines to dynamically change the template matrix of a bundle based off a pre-determined
strategy
"""
class TempStrategy(ABC):

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
    def closing_strat(self, bund):
        pass

class StaticStrat(TempStrategy):

    def open_strat(self, bund):
        return bund

    def close_strat(self,bund):
        return bund
