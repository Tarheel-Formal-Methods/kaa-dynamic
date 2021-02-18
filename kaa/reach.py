import copy
from termcolor import colored

from kaa.timer import Timer
from kaa.bundle import Bundle, BundleTransformer, BundleMode
from kaa.flowpipe import FlowPipe
from kaa.settings import KaaSettings
from kaa.templates import StaticStrat

bolden = lambda string: colored(string, 'white', attrs=['bold'])

class ReachError:

    def __init__(self, msg, total_steps):
        self.message = msg
        self.total_steps = total_steps


"""
Object handling all reachable flowpipe computations.
"""
class ReachSet:

    def __init__(self, model, strat=None, label=""):
        self.model = model
        self.strat = StaticStrat(self.model) if strat is None else strat
        self.flowpipe = FlowPipe(self.model, strat, label)

    """
    Compute reachable set for the alloted number of time steps.
    @params time_steps: number of time steps to carry out the reachable set computation.
            TempStrat: template loading strategy to use during this reachable set computation.
    @returns FlowPipe object containing computed flowpipe
    """
    def computeReachSet(self, time_steps, transmode=BundleMode.AFO):
        transformer = BundleTransformer(self.model, transmode)
        init_box_vol_thres = (100 * self.model.dim) * self.model.bund.getIntersect().volume

        for step in range(time_steps):
            Timer.start('Reachable Set Computation')
            starting_bund = copy.deepcopy(self.flowpipe[step])

            #print("Open: L: {} \n T: {}".format(starting_bund.L, starting_bund.T))
            #print("Open: Offu: {} \n Offl{}".format(starting_bund.offu, starting_bund.offl))

            try:
                self.strat.open_strat(starting_bund, step)
                trans_bund = transformer.transform(starting_bund)
                self.strat.close_strat(trans_bund, step)
            except:
                #if KaaSettings.SaveStateonError:
                #    save_flowpipe_to_disk(self.flowpipe) #Implement this
                #else:
                raise
            
            #print("Close: L: {} \n T: {}".format(trans_bund.L, trans_bund.T))
            #print("Close: Offu: {} Offl{}".format(trans_bund.offu, trans_bund.offl))

            reach_time = Timer.stop('Reachable Set Computation')

            'TODO: Revamp Kaa.log to be output sink handling all output formatting.'
            if not KaaSettings.SuppressOutput:
                print("Computed Step {} -- Time Elapsed: {} sec".format(bolden(step), bolden(reach_time)))

            self.flowpipe.append(trans_bund)

            'Check volume of enveloping box and stop loop if the volume becomes too large.'
            if self.check_reach_size(trans_bund, init_box_vol_thres):
                print("Bundle volume grown to too large of volume. Ending reachable set computation.")
                err = ReachError("Volume too large.", step)
                flowpipe.error = err
                break

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
