import numpy as np

from kaa.templates import TempStrategy
from kaa.experiutil import generate_traj

"""
Local linear approximation strategy.
"""
class LinStrat(TempStrategy):

    def __init__(self, model, iter_steps=2):
        super().__init__(model)
        self.dim = model.dim
        self.iter_steps = iter_steps

        self.unit_dir_mat = np.zeros((self.dim, self.dim))
        for i in range(self.dim):
            self.unit_dir_mat[i][i] = 1

        self.curr_temp = None
        self.curr_dir = None

    def open_strat(self, bund):
        if not self.counter % self.iter_steps:

            #print(f"OPEN DIR/TEMP MAT: {bund.L},  {bund.T}")

            approx_A = self._approx_A(bund, self.dim)
            inv_A = np.linalg.inv(approx_A)
            lin_dir = np.dot(self.unit_dir_mat, inv_A)

            dir_idxs = bund.add_dir(lin_dir)

            if not self.counter:
                bund.add_temp([dir_idxs])
                self.hash_dir('PrevDir', dir_idxs)

            self.unit_dir_mat = lin_dir

        return bund


    def close_strat(self, bund):
        
        if not self.counter % self.iter_steps:
            
            if self.counter:
                bund.remove_dir(self.dir_hash['PrevDir'])
            
        self.counter += 1
        #print(f"CLOSE DIR/TEMP MAT: {bund.L},  {bund.T}")


    def _approx_A(self, bund, num_traj):

        trajs = generate_traj(bund, num_traj, self.iter_steps)
        coeff_mat = np.zeros((self.dim*num_traj,self.dim**2), dtype='float')

        'Initialize the A matrix containing linear constraints.'
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
