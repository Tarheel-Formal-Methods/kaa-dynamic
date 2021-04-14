import numpy as np
from collections import namedtuple
from enum import Enum, auto

from kaa.timer import Timer
from kaa.settings import PlotSettings
from kaa.templates import MultiStrategy
from kaa.linearsystem import calc_box_volume

FlowpipeVolDataTup = namedtuple('FlowpipeVolDataTup', ['FlowpipeConvHullVol', 'FlowpipeEnvelopBoxVol'])
TotalFlowpipeVolTup = namedtuple('FlowpipeVolDataTup', ['TotalFlowpipeConvHullVol', 'TotalFlowpipeEnvelopBoxVol'])

class ReachCompMode(Enum):
    'Computes only the volume statistics of flowpipe. Discards bundle matrices once done.'
    VolMode = auto()

    'Computes only the maximum, minimum values along axes. Discards bundle matrices once done.'
    ProjPlotMode = auto()

    'Most expensive mode in terms of memory usage. Keeps all bundle matrices in memory.'
    PhasePlotMode = auto()


"""
Object encapsulating flowpipe data. A flowpipe in this case will be a sequence of Bundle objects.
"""
class FlowPipe:

    def __init__(self, model, strat, label, reach_comp_mode, flowpipe=None):
        self.flowpipe = [model.bund]
        self.model = model
        self.strat = strat
        self.vars = model.vars
        self.dim = model.dim
        self.label = label
        self.reach_comp_mode = reach_comp_mode


        self.traj_data = None
        self.volume_data = None
        self.proj_data = None
        self.error = None

    """
    Returns a list of strategies which were acting during the reachable set
    computation producing this flowpipe. The strategies are in acting order supplied into
    MultiStrategy.
    """
    @property
    def strats(self):
        return self.strat

    """
    Returns accumlation sum tuple of bundle volumes computed by every volume estimation method.
    """
    @property
    def all_total_volume(self):
        assert self.volume_data, "self.volume_data not initialized. Did you run Flowpipe.calc_flowpipe_data?"
        return TotalFlowpipeVolTup(np.sum(self.volume_data.FlowpipeConvHullVol),
                                   np.sum(self.volume_data.FlowpipeEnvelopBoxVol))

    """
    Return total volume based on most accurate available estimation method.
    Convex hull gets priority for low-dimensioanl systems.
    """
    @property
    def total_volume(self):
        assert self.volume_data, "self.volume_data not initialized. Did you run Flowpipe.calc_flowpipe_data?"
        tot_conv_vol = np.sum(self.volume_data.FlowpipeConvHullVol)

        if tot_conv_vol > 0:
            return tot_conv_vol

        return np.sum(self.volume_data.FlowpipeEnvelopBoxVol)

    @property
    def init_box_volume(self):
        init_box = self.model.init_box
        return calc_box_volume(init_box)

    """
    Wrapper method around self.flowpipe.append function.
    """
    def append(self, bund):
        self.flowpipe.append(bund)

    """
    Crunch requested flowpipe data and delete the bundle list. Memory-conserving procedure.
    """
    def calc_flowpipe_data(self):
        if self.reach_comp_mode == ReachCompMode.VolMode:
            self.volume_data = self.get_volume_data()
            del self.flowpipe

        elif self.reach_comp_mode == ReachCompMode.ProjPlotMode:
            self.proj_data = np.empty((2*self.dim, len(self))

            for var_idx in range(0, 2*self.dim, 2):
                var_min_arr, var_max_arr = self.__calc_bounds_along_var(var_idx)
                self.proj_data[var_idx] = var_min_arr
                self.proj_data[var_idx+1] = var_max_arr

            del self.flowpipe

    def get_strat_flowpipe(self, strat):
        strat_flowpipe = []
        for bund_idx, bund in enumerate(self.flowpipe):
            ptope_strat_list = bund.get_ptopes_by_strat(strat)
            #assert len(ptope_strat_list) != 0, f"Input Strategy must act on bundle object at index {bund_idx}"
            strat_flowpipe.append(ptope_strat_list[0]) #Get first one for now #This stuff should be abstracted and manpulatible.

        return strat_flowpipe

    """
    Returns array of volume data for each bundle in the flowpipe.
    @returns array of volume data.
    """
    def get_volume_data(self, accum=False):
        conv_hull_vol_data = np.empty(len(self.flowpipe))
        envelop_box_vol_data = np.empty(len(self.flowpipe))
        conv_hull_failed = False

        if accum:
            conv_hull_vol_accum = 0
            envelop_box_vol_accum = 0

        for idx, bund in enumerate(self.flowpipe):
            vol_data = bund.getIntersect().volume

            if not vol_data.ConvHullVol:
                conv_hull_failed = True
                conv_hull_vol_data = np.zeros(1) # Empty flowpipe

            if not conv_hull_failed:
                if accum:
                    conv_hull_vol_accum += vol_data.ConvHullVol
                    conv_hull_vol = conv_hull_vol_accum
                else:
                    conv_hull_vol = vol_data.ConvHullVol

                conv_hull_vol_data[idx] = conv_hull_vol

            if accum:
                envelop_box_vol_accum += vol_data.EnvelopBoxVol
                envelop_box_vol = envelop_box_vol_accum
            else:
                envelop_box_vol = vol_data.EnvelopBoxVol

            envelop_box_vol_data[idx] = envelop_box_vol


        return FlowpipeVolDataTup(conv_hull_vol_data, envelop_box_vol_data)


    def __calc_bounds_along_var(self, var_ind):
        Timer.start('Proj')

        'Vector of minimum and maximum points of the polytope represented by parallelotope bundle.'
        y_min, y_max = np.empty(len(self)), np.empty(len(self))

        'Initialize objective function'
        y_obj = np.zeros(self.dim)
        y_obj[var_ind] = 1

        'Calculate the minimum and maximum points through LPs for every iteration of the bundle.'
        for bund_ind, bund in enumerate(self.flowpipe):
            bund_sys = bund.getIntersect()
            y_min[bund_ind] = bund_sys.max_opt(y_obj).fun
            y_max[bund_ind] = bund_sys.min_opt(y_obj).fun
        Timer.stop("Proj")

        return y_min, y_max


    """
    Calculates the flowpipe projection of reachable set against time t.
    @params var: The variable for the reachable set to be projected onto.
    @returns list of minimum and maximum points of projected set at each time step.
    """
    def getProj(self, var_ind):
        assert self.proj_data, "self.proj_data not initialized. Did you run Flowpipe.calc_flowpipe_data?"
        return self.proj_data[2*var_ind], self.proj_data[2*var_ind + 1]

    def __len__(self):
        return len(self.flowpipe)

    def __iter__(self):
        return iter(self.flowpipe)

    def __getitem__(self, index):
        return self.flowpipe[index]

    def __str__(self):
        return self.label
