import numpy as np

from kaa.pykodiak.pykodiak_interface import Kodiak
from kaa.opts.optprod import OptimizationProd

class KodiakProd(OptimizationProd):

    def __init__(self, poly, bund):
        super().__init__(poly, bund)
        self.kodiak_poly = Kodiak.sympy_to_kodiak(self.poly)
    
    def getBounds(self):

        'Unit box bounds'
        bounds = [[0,1] for _ in range(self.bund.sys_dim)]
        jac_mat = np.zeros((self.bund.sys_dim, self.bund.sys_dim))

        lb, ub, _, _ = Kodiak.minmax_diff(self.kodiak_poly, jac_mat, 0, bounds)

        return ub, lb
