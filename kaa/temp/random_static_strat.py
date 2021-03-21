import numpy as np
import random
from math import sin, cos, radians

from kaa.templates import TempStrategy

"""
Strategy which adds random directions and templates at the start of the
reachable set computation.
ONLY WORKS FOR 2D SYSTEMS FOR NOW I.E VANDERPOL.
"""
class RandomStaticStrat(TempStrategy):

    def __init__(self, model, num_ran_temps):
        super().__init__(model)
        self.num_ran_temps = num_ran_temps


    def open_strat(self, bund, step_num):
        if not step_num:
            #ran_dirs_mat = np.empty((self.num_ran_temps, self.dim))

            for row_idx in range(self.num_ran_temps):
                rand_ang_1, rand_ang_2 = [radians(random.randrange(360)), radians(random.randrange(360))]
                ran_dirs = [[cos(rand_ang_1), sin(rand_ang_1)], [cos(rand_ang_2), sin(rand_ang_2)]]

                self.add_ptope_to_bund(bund, ran_dirs, [f"RandDir1Ptope{row_idx}", f"RandDir2Ptope{row_idx}"])


    def close_strat(self, bund, step_num):
        pass

    def reset(self):
        pass
