from kaa.reach import ReachSet
from kaa.plotutil import Plot
from kaa.trajectory import Traj
from kaa.experiutils import get_init_box_borders

class Experiment:

    def __init__(self, model, strat):
        self.model = model
        self.strat = strat

        self.mod_reach = ReachSet(model)
        self.plot = Plot()


    def execute(self, num_steps):

        mod_flow = self.mod_reach.computeReachSet(num_steps, self.strat)
        self.plot.add(mod_flow)

        border_sim_trajs = self.__simulate_border_points()

        for sim_traj in border_sim_trajs:
            self.plot.add(sim_traj)

    def plot_results(self, *var_tup):
        self.plot.plot(*var_tup)

    def __get_init_box(self):

        init_offu = model.bund.offu
        init_offl = model.bund.offl

        return [ [-lower_off, upper_off] for lower_off, upper_off in zip(init_offl, init_offu) ]

    def __simulate_border_points(self, num_steps):
        init_box_inter = self.__get_init_box()
        border_points = get_init_box_borders(init_box_inter)

        trajs = [ Traj(self.model, point, num_steps) for point in border_points ]
        return trajs


class PhaseExperiment(Experiment):

    def __init__(self, model, strat):
        super().__init__(model, strat)

    def plot_results(self, *var_tup):
        self.plot.plot2DPhase(*var_tup)
