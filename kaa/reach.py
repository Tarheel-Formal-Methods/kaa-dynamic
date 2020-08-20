from kaa.timer import Timer
from kaa.bundle import Bundle, BundleTransformer
from kaa.flowpipe import FlowPipe
from kaa.settings import KaaSettings

DefaultStrat = KaaSettings.TempStrat

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
    @returns FlowPipe object with computed flowpipe
    """
    def computeReachSet(self, time_steps, TempStrat=DefaultStrat):

        initial_set = self.model.bund
        transformer = BundleTransformer(self.model.f)
        strat = TempStrat(self.model)
        flowpipe = [initial_set]


        for ind in range(time_steps):
            
            Timer.start('Reachable Set Computation')
            starting_bund = flowpipe[ind]
            opening_bund = strat.open_strat(starting_bund)

            trans_bund = transformer.transform(opening_bund)
            final_bund = strat.close_strat(trans_bund)
            reach_time = Timer.stop('Reachable Set Computation')

            print("Computed Step {0} -- Time Elapsed: {1} sec".format(ind, reach_time))
            flowpipe.append(final_bund)

        return FlowPipe(flowpipe, self.model)
