import matplotlib.pyplot as plt
import numpy as np

from kaa.timer import Timer
from kaa.settings import PlotSettings

import os
"""
Object encapsulating flowpipe data. A flowpipe in this case will be a sequence of bundles.i
"""
class FlowPipe:

    def __init__(self, flowpipe, model, strat):

        self.flowpipe = flowpipe
        self.model = model
        self.strat = strat
        self.vars = model.vars
        self.dim = model.dim

    """
    Calculates the flowpipe projection of reachable set against time t.
    
    @params var: The variable for the reachable set to be projected onto.
    @returns list of minimum and maximum points of projected set at each time step.
    """
    def get2DProj(self, var_ind):
        pipe_len = len(self.flowpipe)

        Timer.start('Proj')
        curr_var = self.vars[var_ind]

        'Vector of minimum and maximum points of the polytope represented by parallelotope bundle.'
        y_min, y_max = np.empty(pipe_len), np.empty(pipe_len)

        'Initialize objective function'
        y_obj = [0 for _ in self.vars]
        y_obj[var_ind] = 1

        'Calculate the minimum and maximum points through LPs for every iteration of the bundle.'
        for bund_ind, bund in enumerate(self.flowpipe):

            bund_sys = bund.getIntersect()

            y_min[bund_ind] = bund_sys.max_opt(y_obj).fun
            y_max[bund_ind] = bund_sys.min_opt(y_obj).fun

        Timer.stop("Proj")

        return y_min, y_max

    @property
    def model_name(self):
        return self.model.name

    def __len__(self):
        return len(self.flowpipe)

    def __iter__(self):
        return iter(self.flowpipe)

    def __str__(self):
        return "{} Len: {}".format(self.strat, len(self))
