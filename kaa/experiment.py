from kaa.reach import ReachSet
from kaa.plotutil import Plot
from kaa.trajectory import Traj
from kaa.experiutil import get_init_box_borders

class Experiment:

    def __init__(self, model, strat=None):
        self.model = model
        self.strat = strat

        self.mod_reach = ReachSet(model)
        self.plot = Plot()


    def execute(self, num_steps):

        border_sim_trajs = self.__simulate_border_points(num_steps)

        for sim_traj in border_sim_trajs:
            self.plot.add(sim_traj)

        mod_flow = self.mod_reach.computeReachSet(num_steps, tempstrat=self.strat)
        self.plot.add(mod_flow)

    def plot_results(self, *var_tup):
        self.plot.plot(*var_tup)

    def __get_init_box(self):

        init_offu = self.model.bund.offu[:self.model.dim] #Assume first dim # of offsets are associated to initial box
        init_offl = self.model.bund.offl[:self.model.dim]

        return [ [-lower_off, upper_off] for lower_off, upper_off in zip(init_offl, init_offu) ]

    def __simulate_border_points(self, num_steps):
        
        init_box_inter = self.__get_init_box()
        border_points = get_init_box_borders(init_box_inter)
        print(border_points)
        trajs = [ Traj(self.model, point, num_steps) for point in border_points ]
        return trajs


class PhasePlotExperiment(Experiment):

    def __init__(self, model, strat=None):
        super().__init__(model, strat)

    def plot_results(self, *var_tup):
        self.plot.plot2DPhase(*var_tup)
