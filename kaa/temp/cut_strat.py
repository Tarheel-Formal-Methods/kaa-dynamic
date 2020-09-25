from kaa.templates import TempStrategy

class CutStrat(TempStrategy):

    def __init__(self, model):
        super().__init__(model)
        self.dim = len(self.model.vars)
        self.counter = 0

    def open_strat(self, bund):

        if self.counter == 100:
            bund.remove_temp([1,2])
            bund.remove_dir([3,4])

        self.counter += 1
        return bund

    def close_strat(self, bund):
        return bund
