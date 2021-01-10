from plotly.offline import plot
from openpyxl import Workbook
import plotly.graph_objects as go
import numpy as np
import os

from kaa.reach import ReachSet
from kaa.plotutil import Plot, TempAnimation
from kaa.trajectory import Traj, TrajCollection
from kaa.settings import PlotSettings, KaaSettings

class Experiment:

    def __init__(self, *inputs, label="Experiment", num_trials=1):
        self.inputs = inputs
        self.plot = Plot()
        self.max_num_steps = 0

        'Assuming that all models use the same dynamics and same initial set for now'
        self.model = inputs[0]['model']
        self.label = label
        self.num_trials = num_trials
        self.results = []

    def execute(self, experi_type="Volume"):
        if experi_type == "Volume":
            for trial_num in range(self.num_trials):
                KaaSettings.RandSeed += trial_num

                flow_labels, flow_vols = zip(*self.gather_vol_data())
                self.results.append(dict(
                    trial_num=trial_num,
                    flow_labels=flow_labels,
                    flow_vols=flow_vols
                ))

    def generate_spreadsheet(self):
        workbook = Workbook()
        sheet = workbook.active

        sheet.append(["Strategy"] + [f"Trial {i+1}" for i in range(self.num_trials)] + ["Mean", "Stdev"])
        'Initialize label-row dictionary'
        row_dict = {label : row_idx + 2 for row_idx, label in enumerate(self.results[0]['flow_labels'])}
        total_num_trials = len(self.results)

        for result in self.results:
            trial_num = result['trial_num']
            for flow_label, flow_vol in zip(result['flow_labels'], result['flow_vols']):
                row = str(row_dict[flow_label])
                sheet['A' + row] = flow_label
                sheet[chr(66 + trial_num) + row] = flow_vol
                sheet[chr(66 + total_num_trials) + row] = f"=AVERAGE(B{row}:{chr(66 + total_num_trials - 1)}{row})"
                sheet[chr(66 + total_num_trials + 1) + row] = f"=STDEV(B{row}:{chr(66 + total_num_trials - 1)}{row})"

        workbook.save(filename=os.path.join(PlotSettings.default_fig_path, self.label + '.xlsx'))
        return workbook

    """
    Execute the reachable set simulations and add the flowpipes to the Plot.
    """
    def gather_plots(self):
        self.output_flowpipes = []
        for experi_input in self.inputs:
            model = experi_input['model']
            strat = experi_input['strat']
            flow_label = experi_input['label']
            num_steps = experi_input['num_steps']

            mod_reach = ReachSet(model, strat=strat, label=flow_label)
            mod_flow = mod_reach.computeReachSet(num_steps)
            #cProfile.runctx('mod_flow = mod_reach.computeReachSet(num_steps, tempstrat=strat, label=flow_label)',None, locals())
            self.plot.add(mod_flow)
            self.output_flowpipes.append(mod_flow)
            self.max_num_steps = max(self.max_num_steps, num_steps)


    def gather_vol_data(self):
        vol_data = []
        for experi_idx, experi_input in enumerate(self.inputs):
            model = experi_input['model']
            strat = experi_input['strat']
            flow_label = experi_input['label']
            num_steps = experi_input['num_steps']

            mod_reach = ReachSet(model, strat=strat, label=flow_label)
            try:
                mod_flow = mod_reach.computeReachSet(num_steps)
                vol_data.append((flow_label, mod_flow.total_volume))
            except Exception as excep:
                vol_data.append((flow_label, str(excep)))

                
        return vol_data

    """
    Plot the results fed into the Plot object
    """
    def plot_results(self, *var_tup, plottrajs=True):
        border_sim_trajs = self.__simulate_border_points(self.max_num_steps)
        if plottrajs:
           self.plot.add(border_sim_trajs)
        self.plot.plot(*var_tup)

    """
    Extract total volume from each experiment given as input.
    """
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
        border_points = __get_init_box_borders(init_box_inter)

        trajs = [Traj(self.model, point, num_steps) for point in border_points]
        return TrajCollection(trajs)

    def __str__(self):
        return self.label

class PhasePlotExperiment(Experiment):

    def __init__(self, *inputs):
        super().__init__(*inputs)

    def plot_results(self, *var_tup, separate=False):
        self.plot.plot2DPhase(*var_tup, separate=separate)

class Animation:

    def __init__(self, experi_input):
        #assert isinstance(experi_input, ExperimentInput), "One ExperimentInput is allowed for animation."
        self.experi_input = experi_input

    def execute(self):
        model = self.experi_input['model']
        strat = self.experi_input['strat']
        label = self.experi_input['label']
        num_steps = self.experi_input['num_steps']

        mod_reach = ReachSet(model, strat=strat)
        mod_flow = mod_reach.computeReachSet(num_steps)
        self.animation = TempAnimation(mod_flow)

    def animate(self, x, y, *strat):
        assert self.animation is not None, "Run Animation.execute first to generate flowpipe to create TempAnimation object."
        self.animation.animate(x, y, *strat)

"""
Wraps around a batch of Experiments for coalesed output retrival.
TODO Get rid of this class. Redundent and can easily be factored into Experiment.
"""
class ExperimentBatch:

    def __init__(self, label=""):
        self.experiments = []
        self.vol_data = []
        self.label = label

    def add_experi(self, experiment):
        assert isinstance(experiment, Experiment), "Only takes Experiment objects."
        self.experiments.append(experiment)

    def execute(self):
        assert len(self.experiments) != 0, "Must add Experiments to ExperimentBatch before executing the batch."
        for experi in self.experiments:
            print(f"\n Executing Experiment {experi.label} \n")
            experi.execute()

    """
    Returns total volume results from consitutent experiments in order which they were added in ExperimentBatch
    """
    def get_vol_data(self):
        return [(str(experi.inputs[0]['strat']), experi.get_total_vol_results()) for experi in self.experiments]

    def get_strat_labels(self):
        return [str(experi.inputs[0]['strat']) for experi in self.experiments]

"""
Find corner vertices for an initial box along with midpoints between the corners.
@params init_box : intervals of the box given as a list of lists where each member's left,right value
                   are the start,end points respectively for the intervals of the box.
@returns list of border points.
"""
def __get_init_box_borders(init_box):

    midpoints = [start + (end - start) / 2 for start, end in init_box]
    border_points = list(product(*init_box))

    for point_idx, point in enumerate(midpoints):
        half_points = [init_inter if point_idx != inter_idx else [point] for inter_idx, init_inter in enumerate(init_box)]
        border_points += list(product(*half_points))

    return border_points

def exec_plot_vol_results(experi, filename):
    labels, vol_data = zip(*experi.gather_vol_only())

    tab_header = dict(values=['Strategy', 'Total Volume'],
                  align='left')
    tab_cells = dict(values=[labels, vol_data],
                  align='left')

    experi_vol_table = go.Table(header=tab_header, cells=tab_cells)

    fig = go.Figure(data=[experi_vol_table])
    fig.write_image(os.path.join(PlotSettings.default_fig_path, filename + 'png'), format='png')
