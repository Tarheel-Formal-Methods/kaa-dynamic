import numpy as np
import random

from kaa.templates import TempStrategy
from kaa.timer import Timer

"""
Strategy which adds random directions and templates at the start of the
reachable set computation.
"""
class RandomStaticStrat(TempStrategy):

    def __init__(self, model, num_ran_temps):
        super().__init__(model)
        self.num_ran_temps = num_ran_temps

    def open_strat(self, bund, step_num):
        if not step_num:
            Timer.start("Random Direction Gen")

            for ptope_idx in range(self.num_ran_temps):
                rand_vec_mat = self.__gen_gaussian_ran_dir()
                rand_vec_labels = [f"RandDir{i}Ptope{ptope_idx}" for i in range(self.dim)]
                self.add_ptope_to_bund(bund, rand_vec_mat, rand_vec_labels)

            Timer.stop("Random Direction Gen")

    def close_strat(self, bund, step_num):
        pass

    def reset(self):
        pass
    """
    Normalize a Gaussian vector to sample surface of n-sphere for higher n
    """
    def __gen_gaussian_ran_dir(self):
        ran_dirs = np.random.normal(size=(self.dim,self.dim))
        norms = np.linalg.norm(ran_dirs, axis=1)

        return ran_dirs / norms

    def __gen_naive_ran_dir(self):
        rand_vec = None
        rand_vec_norm = 2
        while rand_vec_norm > 1:
            rand_vec = [random.uniform(-1, 1) for _ in range(self.dim)]
            rand_vec_norm = sqrt(sum(map(lambda x: x**2, rand_vec)))

        return rand_vec
