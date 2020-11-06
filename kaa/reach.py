import copy
from termcolor import colored 

from kaa.timer import Timer
from kaa.bundle import Bundle, BundleTransformer, BundleMode
from kaa.flowpipe import FlowPipe
from kaa.settings import KaaSettings


DefaultStrat = KaaSettings.DefaultStrat

bolden = lambda string: colored(string, 'white', attrs=['bold'])

"""
Object handling all reachable flowpipe computations.
"""
class ReachSet:

    def __init__(self, model):
        self.model = model

    """
    Compute reachable set for the alloted number of time steps.
    @params time_steps: number of time steps to carry out the reachable set computation.
            TempStrat: template loading strategy to use during this reachable set computation.
    @returns FlowPipe object containing computed flowpipe
    """
    def computeReachSet(self, time_steps, tempstrat=None, transmode=BundleMode.AFO):

        initial_set = self.model.bund
        transformer = BundleTransformer(self.model, transmode)

        strat = tempstrat if tempstrat is not None else DefaultStrat(self.model)
        flowpipe = [initial_set]

        for ind in range(time_steps):
            
            Timer.start('Reachable Set Computation')

            starting_bund = copy.deepcopy(flowpipe[ind])

            #print("Open: L: {} \n T: {}".format(starting_bund.L, starting_bund.T))
            #print("Open: Offu: {} \n Offl{}".format(starting_bund.offu, starting_bund.offl))

            strat.open_strat(starting_bund)
            trans_bund = transformer.transform(starting_bund)
            strat.close_strat(trans_bund)
            
            #print("Close: L: {} \n T: {}".format(trans_bund.L, trans_bund.T))
            #print("Close: Offu: {} Offl{}".format(trans_bund.offu, trans_bund.offl))

            reach_time = Timer.stop('Reachable Set Computation')

            'TODO: Revamp Kaa.log to be output sink handling all output formatting.'
            if not KaaSettings.SuppressOutput:
                print("Computed Step {} -- Time Elapsed: {} sec".format(bolden(ind), bolden(reach_time)))
                
            flowpipe.append(trans_bund)

        return FlowPipe(flowpipe, self.model, strat)
