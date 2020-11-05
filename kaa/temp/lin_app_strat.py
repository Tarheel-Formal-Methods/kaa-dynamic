import numpy as np
from itertools import product
from random import uniform

from kaa.templates import TempStrategy
from kaa.experiutil import generate_traj

"""
Local linear approximation strategy.
"""
class LinStrat(TempStrategy):

    def __init__(self, model, iter_steps=2, cond_threshold=7):
        super().__init__(model)

        self.dim = model.dim
        self.iter_steps = iter_steps

        self.unit_dir_mat = np.zeros((self.dim, self.dim))
        self.__initialize_unit_mat()

        self.cond_threshold = cond_threshold
        self.lin_app_ptope_queue = []

    def open_strat(self, bund):
        if not self.counter % self.iter_steps:

            #print(f"OPEN DIR/TEMP MAT: {bund.L},  {bund.T}"

            approx_A = self.__approx_A(bund, self.dim)
            inv_A = np.linalg.inv(approx_A)
            lin_dir = np.dot(self.unit_dir_mat, inv_A)
            
            #cond_num = np.linalg.cond(lin_dir)
            #print(f"COND NUM: {cond_num}")
            """
            if cond_num > self.cond_threshold:
                norm_lin_dir = self.__normalize_mat(lin_dir)
                print(f"NORM_LIN_DIR: {norm_lin_dir}")

                closest_dirs = self.__find_closest_dirs(norm_lin_dir)
                lin_dir = self.__merge_closest_dirs(norm_lin_dir, closest_dirs)
                print(f"LIN DIR: {lin_dir}")
            """
            new_ptope_labels = [ str((self.counter, dir_idx)) for dir_idx, _ in enumerate(lin_dir) ]
            bund.add_dirs(self, lin_dir,  new_ptope_labels)

            self.lin_app_ptope_queue.append(self.hash_ptope(new_ptope_labels))
            self.unit_dir_mat = lin_dir

        return bund

    def close_strat(self, bund):
        
        if not self.counter % self.iter_steps:
            
            if self.counter:
                self.rm_ptope_from_bund(bund, self.lin_app_ptope_queue.pop(0))

            self.add_ptope_to_bund(bund, self.lin_app_ptope_queue[0])

            
        self.counter += 1
        #print(f"CLOSE DIR/TEMP MAT: {bund.L},  {bund.T}")

    def __approx_A(self, bund, num_traj):

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


    def __merge_closest_dirs(self, dir_mat, closest_dirs):
        print(f"Closest Dirs: {closest_dirs}")
        first_dir, second_dir = (0,1)
        merged_dir = (dir_mat[first_dir] + dir_mat[second_dir]) / 2
        ortho_dir = [ uniform(-1,1) for _ in range(self.dim) ]

        accum_mat = np.delete(dir_mat, [0,1], axis=0)
        accum_mat = np.append(accum_mat, [merged_dir], axis=0)


        'Orthogonalize generated vector through Gram-Schmidt'
        for row in accum_mat:
            ortho_dir -= np.dot(ortho_dir, row) * row

        norm_ortho_dir = ortho_dir / np.linalg.norm(ortho_dir)
        return np.append(accum_mat, [norm_ortho_dir], axis=0)

    def __find_closest_dirs(self, dir_mat):

        closest_pair = None
        closest_dot_prod = 0

        for (first_idx, first_dir), (second_idx, second_dir) in product(enumerate(dir_mat), repeat=2):
            curr_product = np.dot(first_dir, second_dir)

            #print(curr_product)
            if curr_product > closest_dot_prod and first_idx != second_idx:
                closest_pair = (first_idx, second_idx)
                closest_dot_prod = curr_product

        return closest_pair

    def __normalize_mat(self, mat):
        return mat / np.linalg.norm(mat, ord=2, axis=1, keepdims=True)


    def __initialize_unit_mat(self):
    
        for i in range(self.dim):
            self.unit_dir_mat[i][i] = 1

    def __str__(self):
        return "LinApp(Steps:{})".format(self.iter_steps)

"""
    def calc_jacobian(self, point):
        dyns = sp.Matrix(model.f)
        dyns_jac = dyns.jacobian(model.vars)

        jac_func = sp.lambdify(model.vars, dyns_jac, modules='numpy')
        return jac_func(*point)
"""
