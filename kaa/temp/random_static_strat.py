import numpy as np
import random
from math import sin, cos, radians, sqrt

from kaa.templates import TempStrategy

"""
Strategy which adds random directions and templates at the start of the
reachable set computation.
"""
class RandomStaticStrat(TempStrategy):

    def __init__(self, model, num_ran_temps):
        super().__init__(model)
        self.num_ran_temps = num_ran_temps

    """
    Only viable for lower dimensions.
    """
    def open_strat(self, bund, step_num):
        if not step_num:
            for ptope_idx in range(self.num_ran_temps):
                rand_vec_mat = [self.__gen_naive_ran_dir() for _ in range(self.dim)]
                rand_vec_labels = [f"RandDir{i}Ptope{ptope_idx}" for i in range(self.dim)]
                self.add_ptope_to_bund(bund, rand_vec_mat, rand_vec_labels)

    def close_strat(self, bund, step_num):
        pass

    def reset(self):
        pass

    def __gen_naive_ran_dir(self):
        rand_vec = None
        rand_vec_len = 2
        while rand_vec_len > 1:
            rand_vec = [random.uniform(-1, 1) for _ in range(self.dim)]
            rand_vec_len = sqrt(sum(map(lambda x: x**2, rand_vec)))

        return rand_vec
