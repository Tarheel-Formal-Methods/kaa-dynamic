import numpy as np
from collections import namedtuple

from kaa.timer import Timer
from kaa.linearsystem import calc_box_volume

FlowpipeVolDataTup = namedtuple('FlowpipeVolDataTup', ['FlowpipeConvHullVol', 'FlowpipeEnvelopBoxVol'])
TotalFlowpipeVolTup = namedtuple('FlowpipeVolDataTup', ['TotalFlowpipeConvHullVol', 'TotalFlowpipeEnvelopBoxVol'])

"""
Object encapsulating flowpipe data. A flowpipe in this case will be a sequence of Bundle objects.
"""
class FlowPipe:

    def __init__(self, model, strat, label, flowpipe=None):
        self.flowpipe = [model.bund]
        self.num_bunds = 1
        self.model = model
        self.strat = strat
        self.vars = model.vars
        self.dim = model.dim
        self.label = label
        self.error = None
        self.total_comp_time = 0

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
        vol_data = self.get_volume_data()
        return TotalFlowpipeVolTup(np.sum(vol_data.FlowpipeConvHullVol),
                                   np.sum(vol_data.FlowpipeEnvelopBoxVol))

    """
    Return total volume based on most accurate available estimation method.
    Convex hull gets priority for low-dimensioanl systems.
    """
    @property
    def total_volume(self):
        vol_data = self.get_volume_data()
        tot_conv_vol = np.sum(vol_data.FlowpipeConvHullVol)

        if tot_conv_vol > 0:
            return tot_conv_vol

        return np.sum(vol_data.FlowpipeEnvelopBoxVol)

    """
    Returns volume of the initial box defined at the beginning of this flowpipe.
    """
    @property
    def init_box_volume(self):
        init_box = self.model.init_box
        return calc_box_volume(init_box)

    """
    Wrapper method around self.flowpipe.append function.
    """
    def append(self, bund):
        self.flowpipe.append(bund)
        self.num_bunds += 1

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
    def get_volume_data(self, accum=True):
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
        y_min, y_max = np.empty(self.num_bunds), np.empty(self.num_bunds)

        'Initialize objective function'
        y_obj = np.zeros(self.dim)
        y_obj[var_ind] = 1

        'Calculate the minimum and maximum points through LPs for every iteration of the bundle.'
        for bund_ind, bund in enumerate(self.flowpipe):
            bund_sys = bund.getIntersect()
            y_min[bund_ind] = bund_sys.min_obj(y_obj).fun
            y_max[bund_ind] = bund_sys.max_obj(y_obj).fun
        Timer.stop("Proj")

        return y_min, y_max


    """
    Calculates the flowpipe projection of reachable set against time t.
    @params var: The variable for the reachable set to be projected onto.
    @returns list of minimum and maximum points of projected set at each time step.
    """
    def get_proj(self, var_ind):
        return self.__calc_bounds_along_var(var_ind)


    """
    Returns reachable set of the sum of all the variables.
    """
    def get_total_width_reachable_set(self):
        min_width, max_width = np.empty((self.dim, self.num_bunds)), np.empty((self.dim, self.num_bunds))

        for var_idx in range(self.dim):
            var_min, var_max = self.get_proj(var_idx)

            min_width[var_idx] = var_min
            max_width[var_idx] = var_max

        return np.sum(min_width, axis=0), np.sum(max_width, axis=0)

    """
    Returns volume of the box hull of the entire flowpipe.
    """
    def calc_envelop_box_flowpipe_vol(self):
        return calc_box_volume(self.__calc_envelop_box_flowpipe())

    """
    Calculates the box hull of the entire flowpipe.
    """
    def __calc_envelop_box_flowpipe(self):
        flowpipe_envelop_box = [[np.Inf, -np.Inf]] * self.dim
        for bund in self.flowpipe:
            envelop_box = bund.getIntersect().calc_envelop_box()
            for idx, (env_inter, flow_inter) in enumerate(zip(envelop_box, flowpipe_envelop_box)):
                flowpipe_envelop_box[idx][0] = min(env_inter[0], flow_inter[0])
                flowpipe_envelop_box[idx][1] = max(env_inter[1], flow_inter[1])

        return flowpipe_envelop_box

    """
    Calculates the widths of the final bundle in the flowpipe. That is, the width along
    an axis is the maximum value of projection of the bundle along that axis minus its 
    corresponding minimum.
    """
    def calc_final_flowpipe_widths(self):
        final_widths = []
        for var_idx in range(self.dim):
            var_proj_min, var_proj_max = self.get_proj(var_idx)
            final_widths.append(var_proj_max[-1] - var_proj_min[-1])

        return final_widths

    def __len__(self):
        return self.num_bunds

    def __iter__(self):
        return iter(self.flowpipe)

    def __getitem__(self, index):
        return self.flowpipe[index]

    def __str__(self):
        return self.label
