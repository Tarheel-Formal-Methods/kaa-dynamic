import numpy as np
from collections import namedtuple

from kaa.timer import Timer
from kaa.settings import PlotSettings
from kaa.templates import MultiStrategy

FlowpipeVolDataTup = namedtuple('FlowpipeVolDataTup', ['FlowpipeConvHullVol', 'FlowpipeEnvelopBoxVol'])
TotalFlowpipeVolTup = namedtuple('FlowpipeVolDataTup', ['TotalFlowpipeConvHullVol', 'TotalFlowpipeEnvelopBoxVol'])


"""
Object encapsulating flowpipe data. A flowpipe in this case will be a sequence of Bundle objects.
"""
class FlowPipe:

    def __init__(self, model, strat, label, flowpipe=None):
        self.flowpipe = [model.bund] if flowpipe is None else flowpipe
        self.model = model
        self.strat = strat
        self.vars = model.vars
        self.dim = model.dim
        self.label = label
        self.traj_data = None
        self.error = None

    """
    Returns a list of strategies which were acting during the reachable set
    computation producing this flowpipe. The strategies are in acting order supplied into
    MultiStrategy.
    """
    @property
    def strats(self):
        return self.strat.strat_list if isinstance(self.strat, MultiStrategy) else [self.strat]

    @property
    def model_name(self):
        return self.model.name

    """
    Returns accumlation sum of the bundle volumes
    """
    @property
    def total_volume(self):
        vol_tup = self.get_volume_data()

        return TotalFlowpipeVolTup(np.sum(vol_tup.FlowpipeConvHullVol),
                                   np.sum(vol_tup.FlowpipeEnvelopBoxVol))

    """
    Wrapper method around self.flowpipe.append function.
    """
    def append(self, bund):
        self.flowpipe.append(bund)

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

    """
    Calculates the flowpipe projection of reachable set against time t.
    @params var: The variable for the reachable set to be projected onto.
    @returns list of minimum and maximum points of projected set at each time step.
    """
    def getProj(self, var_ind):
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

    def __len__(self):
        return len(self.flowpipe)

    def __iter__(self):
        return iter(self.flowpipe)

    def __getitem__(self, index):
        return self.flowpipe[index]

    def __str__(self):
        return self.label
