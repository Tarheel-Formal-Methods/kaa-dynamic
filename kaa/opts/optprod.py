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
    def getBounds(self):
        pass
