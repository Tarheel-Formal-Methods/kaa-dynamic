import copy

import numpy as np
from tqdm import tqdm

from kaa.timer import Timer
from kaa.bundle import Bundle, BundleTransformer, BundleTransMode
from kaa.flowpipe import FlowPipe
from kaa.settings import KaaSettings
from kaa.templates import StaticStrat, MultiStrategy
from kaa.log import Output


class ReachError:

    def __init__(self, msg, total_steps):
        self.message = msg
        self.total_steps = total_steps


"""
Object handling all reachable flowpipe computations.
"""


class ReachSet:

    def __init__(self, model, strat=None, label="", trans_mode=BundleTransMode.AFO, restrict_inter=None):
        self.model = model
        self.trans_mode = trans_mode
        self.strat = StaticStrat(self.model) if strat is None else strat
        self.flowpipe = FlowPipe(self.model, strat, label)
        self.restrict_inter = restrict_inter if restrict_inter else [-np.Inf, np.Inf]


    """
    Compute reachable set for the alloted number of time steps.
    @params time_steps: number of time steps to carry out the reachable set computation.
            TempStrat: template n  loading strategy to use during this reachable set computation.
    @returns FlowPipe object containing computed flowpipe
    """

    def computeReachSet(self, time_steps):
        transformer = BundleTransformer(self.model, self.trans_mode, self.restrict_inter)
        init_box_vol_thres = ((10 * self.model.dim) * self.model.bund.getIntersect().volume
                              if KaaSettings.UseThreshold
                              else -1)

        with tqdm(range(time_steps)) as step_iter:
            for step in step_iter:
                Timer.start('Reachable Set Computation')
                starting_bund = copy.deepcopy(self.flowpipe[step])  # probably a better way of copying objects.

                try:
                    Timer.start("Open Strategy")
                    self.strat.open_strat(starting_bund, step)
                    Timer.stop("Open Strategy")

                    Timer.start("Bundle Transformation")
                    trans_bund = transformer.transform(starting_bund)
                    Timer.stop("Bundle Transformation")

                    Timer.start("Close Strategy")
                    self.strat.close_strat(trans_bund, step)
                    Timer.stop("Close Strategy")

                except:
                    raise

                reach_time = Timer.stop('Reachable Set Computation')
                self.flowpipe.append(trans_bund)
                step_iter.set_description(f"Step Time: {reach_time}")

                'Check volume of enveloping box and stop loop if the volume becomes too large.'
                if KaaSettings.UseThreshold and self.check_reach_size(trans_bund, init_box_vol_thres):
                    print("Bundle volume grown to too large of volume. Ending reachable set computation.")

                    error = ReachError("Volume too large.", step)
                    self.flowpipe.error = error
                    break

        if not isinstance(self.strat, MultiStrategy):
            self.flowpipe.traj_data = self.strat.fetch_traj_data()

        return self.flowpipe

    def check_reach_size(self, bund, threshold):
        envelop_box_vol = bund.getIntersect().calc_vol_envelop_box()
        return envelop_box_vol > threshold


"""
class FlowpipeSaveLoader:

    @staticmethod
    def save_flowpipe(flowpipe):
        model = flowpipe.model
        with open("")
"""
