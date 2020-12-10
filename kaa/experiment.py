import plotly.graph_objects as go

from kaa.reach import ReachSet
from kaa.plotutil import Plot, TempAnimation
from kaa.trajectory import Traj, TrajCollection
from kaa.experiutil import get_init_box_borders

class ExperimentInput:

    def __init__(self, model, strat=None, label=None):
        self.model = model
        self.strat = strat
        self.label = label

class Experiment:

    def __init__(self, inputs, strat=None):
        self.inputs = inputs if isinstance(inputs, list) else [inputs]
        self.plot = Plot()

        'Assuming that all models use the same dynamics and same initial set for now'
        self.model = inputs[0].model

    """
    Execute the reachable set simulations and add the flowpipes to the Plot.
    """
    def execute(self, num_steps, plottrajs=True):
        border_sim_trajs = self.__simulate_border_points(num_steps)
        self.output_flowpipes = []

        if plottrajs:
            self.plot.add(border_sim_trajs)

        for experi_input in self.inputs:
            model = experi_input.model
            strat = experi_input.strat
            label = experi_input.label

            mod_reach = ReachSet(model)
            mod_flow = mod_reach.computeReachSet(num_steps, tempstrat=strat)
            self.plot.add(mod_flow, label=label)
            self.output_flowpipes.append(mod_flow)

    """
    Plot the results fed into the Plot object
    """
    def plot_results(self, *var_tup):
        self.plot.plot(*var_tup)

    def get_total_vol_results(self):
        assert self.output_flowpipes is not None, "Execute Experiment with ExperimentInputs before retrieving volume data."
        return [flowpipe.total_volume for flowpipe in self.output_flowpipes]

    """
    Extract the initial box intervals from the model
    """
    def __get_init_box(self):
        init_offu = self.model.bund.offu[:self.model.dim] #Assume first dim # of offsets are associated to initial box
        init_offl = self.model.bund.offl[:self.model.dim]

        return [[-lower_off, upper_off] for lower_off, upper_off in zip(init_offl, init_offu)]

    """
    Sample points from the edges of the box and propagate them for a number of steps.
    """
    def __simulate_border_points(self, num_steps):
        init_box_inter = self.__get_init_box()
        border_points = get_init_box_borders(init_box_inter)

        trajs = [Traj(self.model, point, num_steps) for point in border_points]
        return TrajCollection(trajs)

class PhasePlotExperiment(Experiment):

    def __init__(self, inputs):
        super().__init__(inputs)

    def plot_results(self, *var_tup):
        self.plot.plot2DPhase(*var_tup)

class Animation:

    def __init__(self, experi_input):
        assert isinstance(experi_input, ExperimentInput), "One ExperimentInput is allowed for animation."
        self.experi_input = experi_input

    def execute(self, num_steps):
        model = self.experi_input.model
        strat = self.experi_input.strat
        label = self.experi_input.label

        mod_reach = ReachSet(model)
        mod_flow = mod_reach.computeReachSet(num_steps, tempstrat=strat)
        self.animation = TempAnimation(mod_flow)

    def animate(self, x, y, *strat):
        assert self.animation is not None, "Run Animation.execute first to generate flowpipe to create TempAnimation object."
        self.animation.animate(x, y, *strat)

class VolumeExperimentBatch:

    def __init__(self, experiments):
        self.experiments = experiments

    def execute(self, num_steps):
        for experi in self.experiments:
            experi.execute(num_steps)

    def plot_results(self):
        experi_labels = []
        for experi in self.experiments:
            experi_labels.append(f"Box and {str(experi.inputs[0].strat)}") #Only one input to each Experiment. Case checking and sanitation latet.

        #Abstract this out
        experi_vol = [experi.get_total_vol_results() for experi in self.experiments]

        tab_header = dict(values=['Strategy', 'Total Volume'],
                      align='left')
        tab_cells = dict(values=[experi_labels, experi_vol],
                      align='left')

        vol_data = go.Table(header=tab_header, cells=tab_cells)
        fig = go.Figure(vol_data)
        fig.show()
