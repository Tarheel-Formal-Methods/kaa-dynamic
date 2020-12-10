import numpy as np
from itertools import product
from random import uniform

from kaa.templates import TempStrategy

"""
Abstract linear approximation strategy.
"""
class AbstractLinStrat(TempStrategy):

    def __init__(self, model, cond_threshold):
        super().__init__(model)

        self.unit_dir_mat = np.zeros((self.dim, self.dim))
        self.__initialize_unit_mat()

        self.cond_threshold = cond_threshold
        self.lin_app_ptope_queue = []

    def open_strat(self, bund):
        pass

    def close_strat(self, bund):
        pass

    def generate_lin_dir(self, bund):
        approx_A = self.__approx_A(bund, self.dim)
        inv_A = np.linalg.inv(approx_A)
        lin_dir = np.dot(self.unit_dir_mat, inv_A)
        
        cond_num = np.linalg.cond(lin_dir)
        #print(f"COND NUM: {cond_num}")

        if cond_num > self.cond_threshold:
            norm_lin_dir = self.__normalize_mat(lin_dir)
            #print(f"NORM_LIN_DIR: {norm_lin_dir}")

            closest_dirs = self.__find_closest_dirs(norm_lin_dir)
            lin_dir = self.__merge_closest_dirs(norm_lin_dir, closest_dirs)
            #print(f"LIN DIR: {lin_dir}")

        lin_dir_labels = [ str((self.counter, dir_idx)) for dir_idx, _ in enumerate(lin_dir) ]

        return lin_dir, lin_dir_labels

    def __approx_A(self, bund, num_traj):

        trajs = bund.getIntersect().generate_traj(num_traj, self.iter_steps)
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

"""
Local linear approximation strategy.
"""
class LinStrat(AbstractLinStrat):

    def __init__(self, model, iter_steps=2, cond_threshold=7):
        super().__init__(model, cond_threshold)

        self.iter_steps = iter_steps
        self.lin_app_ptope_queue = []

    def open_strat(self, bund):
        if not self.counter % self.iter_steps:
            lin_dir, lin_dir_labels = self.generate_lin_dir(bund)

            ptope_label = self.add_ptope_to_bund(bund, lin_dir, lin_dir_labels)
            self.lin_app_ptope_queue.append(ptope_label)
            self.unit_dir_mat = lin_dir

        return bund

    def close_strat(self, bund):
        if not self.counter % self.iter_steps:
            if self.counter:
                last_ptope =  self.lin_app_ptope_queue.pop(0)
                self.rm_ptope_from_bund(bund, last_ptope)

        self.counter += 1


class DelayedLinStrat(AbstractLinStrat):

    def __init__(self, model, traj_steps=5, num_trajs=100, life_span=1):
        super().__init__(model, traj_steps, num_trajs)

        self.lin_ptope_life = []
        self.life_span = life_span

    def open_strat(self, bund):
        self.__add_new_ptope(bund)

    def close_strat(self, bund):
        'Remove dead templates'
        for ptope_label, life in self.pca_ptope_life:
            if life == 0:
                self.rm_ptope_from_bund(bund, ptope_label)

        self.pca_ptope_life = [(label, life-1) for label, life in self.pca_ptope_life if life > 0]

        #print("After:  L: {} \n T: {}".format(bund.L, bund.T))
        #print("After: offu: {} \n  offl: {}".format(bund.offu, bund.offl))

        self.counter += 1

    def __add_new_ptope(self, bund):
        new_lin_dirs, new_dir_labels = self.generate_lin_dir(bund)
        new_ptope_label = self.add_ptope_to_bund(bund, new_lin_dirs, new_dir_labels)
        self.pca_ptope_life.append((new_ptope_label, self.life_span)) #Add fresh ptope and lifespan to step list

        return new_ptope_label

    def __str__(self):
        return "DelayedPCAStrat-" if self.strat_order is None else f"DelayedPCAStrat{self.strat_order}-"

