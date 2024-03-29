from abc import ABC, abstractmethod

"""
Abstract pass to dictate that every optimization procedure must give an upper and lower bound.
All polynomials passed in will be in sympy's format.
"""
class OptimizationProd(ABC):

    def __init__(self, poly, bund):
        self.poly = poly
        self.bund = bund
        self.vars = bund.vars

    """
    All bounds must be returned as a tuple with the first element being the upper bound and the
    second element being the lower bound.
    """
    @abstractmethod
    def getBounds(self):
        pass

    @abstractmethod
    def getStats(selfs):
        pass

    @abstractmethod
    def __enter__(self):
        return self

    @abstractmethod
    def __exit__(self, exc_type, exc_value, traceback):
        pass
