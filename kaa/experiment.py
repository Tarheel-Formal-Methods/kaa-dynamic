from plotly.offline import iplot
import plotly.graph_objects as go

from kaa.reach import ReachSet
from kaa.plotutil import Plot, TempAnimation
from kaa.trajectory import Traj, TrajCollection
from kaa.experiutil import get_init_box_borders

class Experiment:

    def __init__(self, *inputs, label=""):
        self.inputs = inputs
        self.plot = Plot()
        self.max_num_steps = 0

        'Assuming that all models use the same dynamics and same initial set for now'
        self.model = inputs[0]['model']
        self.label = label

    """
    Execute the reachable set simulations and add the flowpipes to the Plot.
    """
    def execute(self):
        self.output_flowpipes = []

        for experi_input in self.inputs:
            model = experi_input['model']
            strat = experi_input['strat']
            label = experi_input['label']
            num_steps = experi_input['num_steps']

            mod_reach = ReachSet(model)
            mod_flow = mod_reach.computeReachSet(num_steps, tempstrat=strat)
            self.plot.add(mod_flow, label=label)
            self.output_flowpipes.append(mod_flow)
            self.max_num_steps = max(self.max_num_steps, num_steps)

    """
    Plot the results fed into the Plot object
    """
    def plot_results(self, *var_tup, plottrajs=True):
        border_sim_trajs = self.__simulate_border_points(self.max_num_steps)
        if plottrajs:
           self.plot.add(border_sim_trajs)
           
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

    def execute(self):
        model = experi_input['model']
        strat = experi_input['strat']
        label = experi_input['label']
        num_steps = experi_input['num_steps']

        mod_reach = ReachSet(model)
        mod_flow = mod_reach.computeReachSet(num_steps, tempstrat=strat)
        self.animation = TempAnimation(mod_flow)

    def animate(self, x, y, *strat):
        assert self.animation is not None, "Run Animation.execute first to generate flowpipe to create TempAnimation object."
        self.animation.animate(x, y, *strat)

class ExperimentBatch:

    def __init__(self, experiments=[], label=""):
        assert isinstance(experiments, list), "experiments keyword must contain list of Experiment objects."
        self.experiments = experiments
        self.label = label

    def add_experi(self, experiment):
        assert isinstance(experiment, Experiment), "Only takes Experiment objects."
        self.experiments.append(experiment)

    def execute(self):
        assert len(self.experiments) != 0, "Must add Experiments to ExperimentBatch before executing the batch."
        for experi in self.experiments:
            print(f"\n Executing Experiment {experi.label} \n")
            experi.execute()

    def get_vol_data(self):
        return [experi.get_total_vol_results() for experi in self.experiments]

    def get_strat_labels(self):
        return [str(experi.inputs[0]['strat']) for experi in self.experiments]
        
def exec_plot_vol_results(*experi_bat):

    experi_tables = []
    for experi_bat_idx, experi_bat in enumerate(experi_bat):
        experi_bat.execute()
        
        tab_header = dict(values=['Strategy', 'Total Volume'],
                  align='left')
        tab_cells = dict(values=[experi_bat.get_strat_labels(),experi_bat.get_vol_data()],
                  align='left')

        experi_vol_table = go.Table(header=tab_header, cells=tab_cells)
        experi_tables.append(experi_vol_table)

    fig = go.Figure(data=experi_tables)
    iplot(fig)
