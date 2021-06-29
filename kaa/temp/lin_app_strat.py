import numpy as np
from itertools import product
import random as rand

from kaa.templates import TempStrategy, GeneratedDirs
from kaa.settings import KaaSettings
from kaa.log import Output
from kaa.timer import Timer

"""
Abstract linear approximation strategy.
"""


class AbstractLinStrat(TempStrategy):

    def __init__(self, model, num_trajs, cond_threshold):
        super().__init__(model)
        self.accum_lin_trans = initialize_unit_mat(self.dim)
        self.cond_threshold = cond_threshold
        self.num_trajs = num_trajs

    def open_strat(self, bund):
        pass

    def close_strat(self, bund):
        pass

    """
    Generate Linear Approximation directions from Bundle object or fetch them from
    GeneratedLinDirs object.
    @params bund: Bundle object
            step_num: step number of computation
    @returns tuple of matrix containing linear app. directions and
             labels for each of those directions.
    """

    def generate_lin_dir(self, bund, step_num):
        if self.dirs is None:
            Timer.start("Linear Direction Generation")
            inv_A = self.__approx_inv_A(bund)  # Approx inverse linear transform.

            lin_dir_mat = np.matmul(self.accum_lin_trans, inv_A)
            if KaaSettings.NormalizeLinDir:
                cond_num = np.linalg.cond(lin_dir_mat)

                if cond_num > self.cond_threshold:
                    Output.prominent("Using Normalization method.")
                    norm_lin_dir = normalize_mat(lin_dir_mat)
                    closest_dirs = find_closest_dirs(norm_lin_dir)
                    lin_dir_mat = merge_closest_dirs(norm_lin_dir, closest_dirs, self.dim)

            Timer.stop("Linear Direction Generation")
            self.accum_lin_trans = lin_dir_mat  # This makes senese for now as we assume we generate directions once per step.
        else:
            lin_dir_mat = self.dirs.get_dirs_at_step(step_num)

        lin_dir_labels = [str((step_num, dir_idx)) for dir_idx, _ in enumerate(lin_dir_mat)]

        return lin_dir_mat, lin_dir_labels

    """
    Approximates the closest linear transformation from a sample of trajectories starting from the
    within the bundle and ending some steps propagated forward.
    @params bund: Bundle object.
    @returns best-fit linear transformation matrix.
    """

    def __approx_inv_A(self, bund):

        trajs = self.generate_trajs(bund, self.num_trajs)
        start_end_tup = [(t.start_point, t.end_point) for t in trajs]

        try:
            approx_A = approx_lin_trans(start_end_tup, self.dim)
            # print(f"LEAST SQ Matrix Output: {approx_A}")

            inv_A = np.linalg.inv(approx_A)
            # print(f"INV Matrix Output: {inv_A}")

        except np.linalg.LinAlgError:
            Output.warning("USING LEAST SQ INVERSE DUE TO SINGULAR VALUE ERROR")
            inv_A = approx_inv_lin_trans(start_end_tup, self.dim)

        return inv_A


"""
Local linear approximation strategy.
"""


class LinStrat(AbstractLinStrat):

    def __init__(self, model, iter_steps=1, num_trajs=-1, cond_threshold=7):

        num_trajs = 2 * model.dim if num_trajs < 0 else num_trajs
        super().__init__(model, num_trajs, cond_threshold)
        self.iter_steps = iter_steps
        self.ptope_queue = []

    """
    Opening LinApp routine
    """

    def open_strat(self, bund, step_num):
        if (self.iter_steps and
                not step_num % self.iter_steps):
            lin_dir, lin_dir_labels = self.generate_lin_dir(bund, step_num)  # Generate or fetch the LinApp directions.

            'Add the components to the bundle and save last generated ptope data.'
            ptope_label = self.add_ptope_to_bund(bund, lin_dir, lin_dir_labels)
            self.ptope_queue.append(ptope_label)

        return bund

    """
    Closing LinApp routine
    """

    def close_strat(self, bund, step_num):
        if (self.iter_steps > 0 and
                not step_num % self.iter_steps and
                step_num > 0):
            self.rm_ptope_from_bund(bund, self.ptope_queue.pop(0))

    """
    Reset the strategy for a new round of computation.
    """

    def reset(self):
        self.ptope_queue = []
        self.accum_lin_trans = initialize_unit_mat(self.dim)

    def __str__(self):
        return f"LinApp(Steps:{self.iter_steps})"


class SlidingLinStrat(AbstractLinStrat):

    def __init__(self, model, lifespan=1, num_trajs=-1, cond_threshold=7):
        num_trajs = 2 * model.dim if num_trajs < 0 else num_trajs
        super().__init__(model, num_trajs, cond_threshold)
        self.lin_ptope_life_counter = {}
        self.lifespan = lifespan

    """
    Opening LinApp routine
    """

    def open_strat(self, bund, step_num):
        # Output.bold_write(f"Calling Open Strat for {str(self)}")
        if self.lifespan > 0:
            self.__add_new_ptope(bund, step_num)

            'Remove dead templates'
            for ptope_label in list(self.lin_ptope_life_counter.keys()):
                self.lin_ptope_life_counter[ptope_label] -= 1

                if self.lin_ptope_life_counter[ptope_label] == 0:
                    self.rm_ptope_from_bund(bund, ptope_label)
                    self.lin_ptope_life_counter.pop(ptope_label)

            assert len(self.lin_ptope_life_counter) == min(step_num + 1,
                                                           self.lifespan), f"Number of templates don't match. {len(self.lin_ptope_life_counter)} at step num {step_num}"

    """
    Closing LinApp routine
    """

    def close_strat(self, bund, step_num):
        # Output.bold_write(f"Calling Closing Strat for {str(self)}")
        pass

    """
    Auxiliary method to generate LinApp. directions and add them as
    directions for a new ptope every step.
    @params bund: Bundle object
            step_num: step number in computation
    @returns: label string for generated ptope
    """

    def __add_new_ptope(self, bund, step_num):
        new_lin_dirs, new_dir_labels = self.generate_lin_dir(bund, step_num)
        new_ptope_label = self.add_ptope_to_bund(bund, new_lin_dirs, new_dir_labels)
        self.lin_ptope_life_counter[new_ptope_label] = self.lifespan + 1  # Add fresh ptope and lifespan to step list

        return new_ptope_label

    """
    Reset the strategy for a new round of computation.
    """

    def reset(self):
        self.lin_ptope_life_counter = {}
        self.accum_lin_trans = initialize_unit_mat(self.dim)

    def __str__(self):
        return f"SlidingLinStrat(Steps:{self.lifespan})"


class GeneratedLinDirs(GeneratedDirs):

    def __init__(self, model, num_steps, num_trajs, cond_threshold=7, dir_mat=None, sampled_points=None):

        if dir_mat is None or sampled_points is None:  # dir_mat is set if pre-generated directions are used during computation.
            self.unit_dir_mat = initialize_unit_mat(model.dim)
            self.cond_threshold = cond_threshold
            self.num_trajs = num_trajs
            lin_dirs, traj_mat = self.__generate_lin_dir(model, num_steps)

            super().__init__(model, lin_dirs, traj_mat)
        else:
            super().__init__(model, dir_mat, sampled_points)

    """
    Generates the linear approximation directions based on trajectories taken starting from
    initial box.
    @params: model: Model object
             num_steps: number of steps to pre-generate directions. These directions
                        will be used in the reachable set computation involving the input model
                        object.
    @returns matrix of generated directions
    """

    def __generate_lin_dir(self, model, num_steps):
        bund = model.bund
        dim = model.dim

        generated_lin_dir_mat = np.empty((dim * num_steps, dim))

        trajs = bund.getIntersect().generate_ran_trajs(self.num_trajs, num_steps)  # trajs is TrajCollecton object'

        for step in range(num_steps):
            start_end_tup = [(t[step], t[step + 1]) for t in trajs]  # Get all trajectory points at step t

            try:
                approx_A = approx_lin_trans(start_end_tup, dim)
                inv_A = np.linalg.inv(approx_A)

            except np.linalg.LinAlgError:
                Output.warning("USING LEAST SQ INVERSE DUE TO SINGULAR VALUE ERROR")
                inv_A = approx_inv_lin_trans(start_end_tup, dim)

            lin_dir = np.matmul(self.unit_dir_mat, inv_A)
            cond_num = np.linalg.cond(lin_dir)

            if cond_num > self.cond_threshold:
                norm_lin_dir = normalize_mat(lin_dir)
                closest_dirs = find_closest_dirs(norm_lin_dir)
                lin_dir = merge_closest_dirs(norm_lin_dir, closest_dirs, dim)

            generated_lin_dir_mat[step * dim: (step + 1) * dim, ] = lin_dir
            self.unit_dir_mat = lin_dir

        return generated_lin_dir_mat, trajs.get_mat()


"""
Approximate the linear transformation.
"""


def approx_lin_trans(dom_ran_tup, dim):
    return __least_sq_trans(dom_ran_tup, dim)


"""
Approximate the inverse linear transformation.
"""


def approx_inv_lin_trans(dom_ran_tup, dim):
    swapped_dom_ran_tup = [(ran, dom) for dom, ran in dom_ran_tup]
    return __least_sq_trans(swapped_dom_ran_tup, dim)


"""
Routine to approximate the best-fit linear transformation matrix which matches the data given in domain-range tuple argument.
Uses np.linalg.lstsq to accomplish this.
@params: dom_ran_tup: tuple of points in system such that the first index is the point x
                      and second index is Ax for some desired linear transformation A
         dim: dimension of system
@returns best-fit linear transformation matrix.
"""


def __least_sq_trans(dom_ran_tup, dim):
    num_data_points = len(dom_ran_tup)
    coeff_mat = np.zeros((dim * num_data_points, dim ** 2), dtype='float')

    'Initialize the A matrix containing linear constraints.'
    for t_idx, t in enumerate(dom_ran_tup):
        end_point = t[0]
        for i in range(dim):
            for j in range(dim):
                coeff_mat[i + dim * t_idx][i * dim + j] = end_point[j]

    b_mat_vec = np.asarray([t[1] for t in dom_ran_tup], dtype='float')
    b_mat = b_mat_vec.flatten()

    m = np.linalg.lstsq(coeff_mat, b_mat, rcond=None)[0]
    return m.reshape((dim, dim))


"""
Takes the directions which are most similar to each other, as determined by the inner product, and
merges the two by taking the average of the two vectors. The additional missing vector is replaced by one
which is orthogonal to resulting set of vectors.
@params dir_mat: matrix of direction vectors as rows.
        closest_dirs: tuple of row indices in dir_mat indictaing two closest direction vectors
        dim: dimension of system.
@returns new directions matrix with merged directions
"""


def merge_closest_dirs(dir_mat, closest_dirs, dim):
    randgen = rand.Random(KaaSettings.RandSeed)
    first_dir, second_dir = closest_dirs
    merged_dir = (dir_mat[first_dir] + dir_mat[second_dir]) / 2
    ortho_dir = [randgen.uniform(-1, 1) for _ in range(dim)]

    accum_mat = np.delete(dir_mat, [0, 1], axis=0)
    accum_mat = np.append(accum_mat, [merged_dir], axis=0)

    'Orthogonalize generated vector through Gram-Schmidt'
    for row in accum_mat:
        ortho_dir -= np.dot(ortho_dir, row) * row

    norm_ortho_dir = ortho_dir / np.linalg.norm(ortho_dir)
    return np.append(accum_mat, [norm_ortho_dir], axis=0)


"""
Finds the two closest direction vectors by taking pairwise inner product between all vectors.
Returns the indices of the two closest ones.
@params dir_mat: directions matrix
@returns tuple of indices of closest pair
"""


def find_closest_dirs(dir_mat):
    closest_pair = None
    closest_dot_prod = 0

    for (first_idx, first_dir), (second_idx, second_dir) in product(enumerate(dir_mat), repeat=2):
        curr_product = np.dot(first_dir, second_dir)

        if abs(curr_product) > closest_dot_prod and first_idx != second_idx:
            closest_pair = (first_idx, second_idx)
            closest_dot_prod = curr_product

    return closest_pair


"""
Normalizes the row of input matrix
@params mat: matrix
@returns normalized matrix.
"""


def normalize_mat(mat):
    return mat / np.linalg.norm(mat, ord=2, axis=1, keepdims=True)


def initialize_unit_mat(dim):
    mat = np.zeros((dim, dim))
    for i in range(dim):
        mat[i][i] = 1
    return mat
