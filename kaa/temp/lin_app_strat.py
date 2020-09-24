import sympy as sp
import numpy as np

from kaa.templates import TempStrat

class LinStrat(TempStrat):

    def __init__(self, model):
        super().__init__(model)
        self.counter = 0

    def open_strat(self, bund):
        if not self.counter % 2:

            if self.counter:
                bund.remove_temp(-1)
                bund.remove_dir(self.dir_hash['invtemp'])

            local_jac = self.calc_jacobian(point)
            inv_jac = np.linalg.inv(local_jac)

            bund.add_dir(inv_jac)


    def close_strat(self):
        pass

    def calc_jacobian(self, point):
        dyns = sp.Matrix(model.f)
        dyns_jac = dyns.jacobian(model.vars)

        jac_func = sp.lambdify(model.vars, dyns_jac, modules='numpy')
        return jac_func(*point)
