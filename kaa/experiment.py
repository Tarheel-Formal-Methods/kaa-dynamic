from plotly.offline import plot
from openpyxl import Workbook
import plotly.graph_objects as go
import numpy as np
import os

from kaa.reach import ReachSet
from kaa.plotutil import Plot, TempAnimation
from kaa.trajectory import Traj, TrajCollection
from kaa.settings import PlotSettings, KaaSettings
from kaa.templates import MultiStrategy
from kaa.temp.pca_strat import AbstractPCAStrat, GeneratedPCADirs
from kaa.temp.lin_app_strat import AbstractLinStrat, GeneratedLinDirs
from kaa.timer import Timer

class SpreadSheet:

    def __init__(self, workbook, row_dict):
        self.workbook = workbook
        self.row_dict = row_dict

class Experiment:

    def __init__(self, *inputs, label="Experiment", num_trials=1):
        self.inputs = inputs
        self.plot = Plot()
        self.max_num_steps = 0

        'Assuming that all models use the same dynamics and same initial set for now'
        self.model = inputs[0]['model']
        self.label = label
        self.num_trials = num_trials

    """
    Execute experiment and dump results into spreadsheet.
    """
    def execute(self, experi_type="Volume"):
        if experi_type == "Volume":
            spreadsheet = self.__generate_sheet()

            for experi_input in self.inputs:
                experi_strat = experi_input['strat']
                experi_num_steps = experi_input['num_steps']
                experi_num_trajs = experi_input['num_trajs']

                self.__generate_dirs(experi_strat, experi_num_steps, experi_num_trajs)

                for trial_num in range(self.num_trials):
                    print(f"RUNNING EXPERIMENT {experi_input['label']} TRIAL:{trial_num}")
                    flow_label, flow_vol = self.__gather_vol_data(experi_input)
                    self.__save_data_into_sheet(spreadsheet, trial_num, flow_label, flow_vol)

                    self.__update_seed()
                    experi_strat.reset()
                    self.__generate_dirs(experi_strat, experi_num_steps, experi_num_trajs)

    """
    Saves data into a desired cell in spreadsheet.
    """
    def __save_data_into_sheet(self, spreadsheet, trial_num, flow_label, data):
        workbook = spreadsheet.workbook
        row_dict = spreadsheet.row_dict

        column_offset = trial_num
        row_offset = row_dict[flow_label]

        sheet = workbook.active
        sheet[chr(66 + column_offset) + str(row_offset)] = data

        if column_offset == self.num_trials - 1:
            sheet[chr(66 + self.num_trials) + str(row_offset)] = f"=AVERAGE(B{row_offset}:{chr(66 + self.num_trials - 1)}{row_offset})"
            sheet[chr(66 + self.num_trials + 1) + str(row_offset)] = f"=STDEV(B{row_offset}:{chr(66 + self.num_trials - 1)}{row_offset})"

        workbook.save(filename=os.path.join(PlotSettings.default_fig_path, self.label + '.xlsx'))

    """
    Pre-generates directions based on the current global seed for Random.
    """
    def __generate_dirs(self, strat, num_steps, num_trajs):
        if isinstance(strat, MultiStrategy):
            for st in strat.strats:
                self.__generate_dirs_by_strat(st, num_steps, num_trajs)
        else:
            self.__generate_dirs_by_strat(strat, num_steps, num_trajs)

    """
    Auxiliary method to generate directions based on strategy type.
    """
    def __generate_dirs_by_strat(self, strat, num_steps, num_trajs):
        if isinstance(strat, AbstractPCAStrat):
            strat.pca_dirs = GeneratedPCADirs(strat.model, num_steps, num_trajs)
        elif isinstance(strat, AbstractLinStrat):
            strat.lin_dirs = GeneratedLinDirs(strat.model, num_steps, num_trajs)
        else:
            raise RuntimeError("Strategies have to be a PCAStrat or a LinStrat")

    """
    Initializes openpyxl spreadsheet to dump resulting data.
    """
    def __generate_sheet(self):
        workbook = Workbook()
        sheet = workbook.active
        sheet.append(["Strategy"] + [f"Trial {i+1}" for i in range(self.num_trials)] + ["Mean", "Stdev"])

        'Initialize label-row dictionary'
        row_dict = {experi_input['label'] : row_idx + 2 for row_idx, experi_input in enumerate(self.inputs)}

        for experi_input in self.inputs:
            flow_label = experi_input['label']
            row = row_dict[flow_label]
            sheet['A' + str(row)] = flow_label

        workbook.save(filename=os.path.join(PlotSettings.default_fig_path, self.label + '.xlsx'))
        return SpreadSheet(workbook, row_dict)

    """
    Update global random seed.
    """
    def __update_seed(self, offset=1):
        KaaSettings.RandSeed += offset

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


    def __gather_vol_data(self, experi_input):
        model = experi_input['model']
        strat = experi_input['strat']
        flow_label = experi_input['label']
        num_steps = experi_input['num_steps']

        mod_reach = ReachSet(model, strat=strat, label=flow_label)
        try:
            mod_flow = mod_reach.computeReachSet(num_steps)
            return (flow_label, mod_flow.total_volume)
        except Exception as excep:
            raise
            return (flow_label, str(excep))

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
