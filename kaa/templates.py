import numpy as np
from abc import ABC, abstractmethod

"""
Object containing routines to dynamically change the template matrix of a bundle based off a pre-determined
strategy
"""
class TempStrategy(ABC):

    def __init__(self, model):
        self.model = model

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
