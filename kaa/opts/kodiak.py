import numpy as np

from kaa.pykodiak.pykodiak_interface import Kodiak
from kaa.opts.optprod import OptimizationProd
from kaa.log import Output

class KodiakProd(OptimizationProd):

    def __init__(self, poly, bund):
        super().__init__(poly, bund)
        self.kodiak = Kodiak()

        for var in bund.vars:
            self.kodiak.add_variable(str(var))
        self.kodiak_poly = self.kodiak.sympy_to_kodiak(self.poly)

    def getBounds(self):
        'Unit box bounds'
        bounds = [[0,1] for _ in range(self.bund.dim)]
        jac_mat = np.zeros((self.bund.dim, self.bund.dim))

        #Output.bold_write("Calling Kodiak")
        #Output.write(f"INPUT POLY: {self.poly}")
        lb, ub, _, _ = self.kodiak.minmax_diff(self.kodiak_poly, jac_mat, 0, bounds)
        #Output.bold_write("Out of Kodiak")

        return ub, lb
