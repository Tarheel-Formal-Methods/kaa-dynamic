import sympy as sp
import numpy as np

from functools import reduce
from math import comb
from operator import mul,add
from itertools import product
from kaa.opts.optprod import OptimizationProd

class BernsteinProd(OptimizationProd):

    def __init__(self, poly, bund):
        super().__init__(poly, bund)
        self.poly = sp.Poly(poly, self.vars)
        self.var_num = len(self.vars)
        self.degree = self._getDegree()

        # List of all relevant monomials required to calculate the Bernstein coefficients.
        self.monom_list = self._computeLowerDeg(self.degree)
        self.bern_coeff = np.empty(len(self.monom_list))

        for idx, monom_deg in enumerate(self.monom_list):
            self.bern_coeff[idx] = self._computeIthBernCoeff(monom_deg)

    """
    Computes and returns the maximum and minimum Bernstein coefficients for self.poly.
    """
    def getBounds(self):
        return max(self.bern_coeff), min(self.bern_coeff)

    def getStats(self):
        return self.getBoundIndices()

    """
    Determines indices of monomials where maximum/minimum coefficients are locateed.
    """
    def getBoundIndices(self):
        max_idxs = np.argmax(self.bern_coeff, axis=0)
        min_idxs = np.argmin(self.bern_coeff, axis=0)

        if isinstance(max_idxs, np.int64):
            max_idxs = [max_idxs]
        if isinstance(min_idxs, np.int64):
            min_idxs = [min_idxs]

        max_monom_deg = [str(self.monom_list[max_idx]) for max_idx in max_idxs]
        min_monom_deg = [str(self.monom_list[min_idx]) for min_idx in min_idxs]

        return max_monom_deg, \
               min_monom_deg
    """
    Computes and returns the ith Bernstein coefficient.
    @params i: the degree of the desired Bernstein polynomial.
    """
    def _computeIthBernCoeff(self, i):
        lower_degs = self._computeLowerDeg(i)

        bern_sum_list = []
        for j in lower_degs:
            coeff_mul_list = []

            for ind, _ in enumerate(self.degree):
                j_coeff = comb(i[ind], j[ind]) / comb(self.degree[ind], j[ind])
                coeff_mul_list.append(j_coeff)

            bern_coef = reduce(mul, coeff_mul_list)
            poly_coef = self.poly.coeff_monomial(self._getMonom(j))
            bern_sum_list.append(bern_coef * poly_coef)

        return reduce(add, bern_sum_list)

    """
    Computes the degrees lower than supplied index.
    Returns a list of all lower degrees.
    @params i: desired degree
    """
    def _computeLowerDeg(self, i):
        assert len(i) == len(self.degree)

        iterators = [range(idx+1) for idx in i]
        return list(product(*iterators))

    """
    Returns the total degree of self.poly
    """
    def _getDegree(self):
        monom_tups = self.poly.monoms()
        degree = []

        for var_index, _ in enumerate(self.vars):
             var_deg = max([monom[var_index] for monom in monom_tups])
             degree.append(var_deg)

        return degree

    """
    Returns the sympy monomial of specific degree.
    @params j: degree of monomial.
    """
    def _getMonom(self, j):
        var_monom = [ self.vars[i]**j[i] for i in range(self.var_num) ]
        expr = reduce(mul, var_monom)

        monomial = sp.Poly(expr, self.vars)
        return monomial.as_expr()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        pass
