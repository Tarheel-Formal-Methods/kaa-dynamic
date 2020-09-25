import sympy as sp
import numpy as np

from kaa.templates import TempStrategy
from kaa.experiutil import generate_traj

class LinStrat(TempStrategy):

    def __init__(self, model, iter_steps=2):
        super().__init__(model)
        self.iter_steps = iter_steps
        self.counter = 0
        
    def open_strat(self, bund):
        if not self.counter % self.iter_steps:

            if self.counter:
                bund.remove_temp(-1)
                bund.remove_dir(self.dir_hash['LinDir'])
                self.pop_hash_dir('LinDir')

            approx_A = self._approx_A(bund, self.dim)
            inv_A = np.linalg.inv(approx_A)

            unit_dir_mat = bund.L[:self.dim]
            lin_dir = np.dot(unit_dir_mat, inv_A)

            dir_idxs = bund.add_dir(lin_dir)
            bund.add_temp([dir_idxs])
            self.hash_dir('LinDir', dir_idxs)

            print("Calculated")

        self.counter += 1
        # print("Calculated")
        return bund


    def close_strat(self, bund):
        pass

    def _approx_A(self, bund, num_traj):
        trajs = generate_traj(bund, num_traj, self.iter_steps)
        coeff_mat = np.zeros((self.dim*num_traj,self.dim**2), dtype='float')
        
        for t_idx, t in enumerate(trajs):
            for i in range(self.dim):
                for j in range(self.dim):
                    coeff_mat[i+self.dim*t_idx][i*self.dim+j] = t.start_point[j]

        b_mat_vec = np.asarray([t.end_point for t in trajs], dtype='float')
        b_mat = b_mat_vec.flatten()

        m = np.linalg.lstsq(coeff_mat, b_mat, rcond=None)[0]
        return m.reshape((self.dim,self.dim))

    def __str__(self):
        return "LinApp(Steps:{})".format(self.iter_steps)

"""
    def calc_jacobian(self, point):
        dyns = sp.Matrix(model.f)
        dyns_jac = dyns.jacobian(model.vars)

        jac_func = sp.lambdify(model.vars, dyns_jac, modules='numpy')
        return jac_func(*point)
"""
