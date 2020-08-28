import numpy as np

from kaa.templates import TempStrategy, Template
from kaa.bundle import Bundle


class SimpleStrat(TempStrategy):

    def __init__(self, model):
        super().__init__(model)
        self.toggle = True
        self.last_temp = np.copy(model.bund.T[1:])
        self.iter = 0
        self.change_rate = 3

        #print("last temp: {}".format(self.last_temp))

    def open_strat(self, bund):
        T = Template(bund)
        if self.toggle:
            T.remove(-1)
            T.remove(-1)
        else:
            T.add(self.last_temp)

        self.toggle = not self.toggle
        return Bundle(T.M, bund.L, bund.offu, bund.offl, bund.vars)
        
    def close_strat(self, bund):
        return bund
