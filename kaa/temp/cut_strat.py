from kaa.templates import TempStrategy


"""
Strategy which removes the non-trivial templates from SIR specifically.
Mostly created as a sanity test to visualize how dynamic template modification would affect
the reachable set.
"""
class CutStrat(TempStrategy):

    def __init__(self, model):
        super().__init__(model)
        self.dim = len(self.model.vars)
        self.counter = 0

    def open_strat(self, bund):

        'Remove the non-trivial templates at time step 100'
        if self.counter == 100:
            bund.remove_temp([1,2])
            bund.remove_dir([3,4])

        self.counter += 1
        return bund

    def close_strat(self, bund):
        return bund
