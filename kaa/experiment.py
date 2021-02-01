from plotly.offline import plot
from openpyxl import Workbook
from abc import ABC, abstractmethod
import plotly.graph_objects as go
import numpy as np
import os

from kaa.reach import ReachSet
from kaa.plotutil import Plot, TempAnimation
from kaa.trajectory import Traj, TrajCollection
from kaa.settings import PlotSettings, KaaSettings
from kaa.templates import MultiStrategy, GeneratedDirs
from kaa.temp.pca_strat import AbstractPCAStrat, GeneratedPCADirs
from kaa.temp.lin_app_strat import AbstractLinStrat, GeneratedLinDirs
from kaa.log import Output
from kaa.timer import Timer

class SpreadSheet:

    def __init__(self, workbook, row_dict):
        self.workbook = workbook
        self.row_dict = row_dict

class DirSaveLoader:

    @staticmethod
    def load_dirs(model, num_steps, num_trajs, seed):
        pca_dirs_from_file = DirSaveLoader.load_and_reshape(model, os.path.join(KaaSettings.DataDir, f"PCA{model}(T:{num_trajs})(Steps:{num_steps})(Seed:{seed}).txt"))
        lin_dirs_from_file = DirSaveLoader.load_and_reshape(model, os.path.join(KaaSettings.DataDir, f"Lin{model}(T:{num_trajs})(Steps:{num_steps})(Seed:{seed}).txt"))

        pca_gendir_obj_list = DirSaveLoader.wrap_pca_dirs(model, pca_dirs_from_file)
        lin_gendir_obj_list = DirSaveLoader.wrap_lin_dirs(model, lin_dirs_from_file)

        return list(zip(pca_gendir_obj_list, lin_gendir_obj_list))

    @staticmethod
    def save_dirs(model, num_steps, num_trajs, seed, gen_dirs_list):
        pca_dir_stack, lin_dir_stack = zip(*gen_dirs_list)

        pca_dir_stack = [gen_dirs.dir_mat for gen_dirs in pca_dir_stack]
        lin_dir_stack = [gen_dirs.dir_mat for gen_dirs in lin_dir_stack]

        pca_data_to_save = np.column_stack([dirs.flatten() for dirs in pca_dir_stack])
        lin_data_to_save = np.column_stack([dirs.flatten() for dirs in lin_dir_stack])

        np.savetxt(os.path.join(KaaSettings.DataDir, f"PCA{model}(T:{num_trajs})(Steps:{num_steps})(Seed:{seed}).txt"), pca_data_to_save, delimiter=',')
        np.savetxt(os.path.join(KaaSettings.DataDir, f"Lin{model}(T:{num_trajs})(Steps:{num_steps})(Seed:{seed}).txt"), lin_data_to_save, delimiter=',')

    def load_and_reshape(model, datapath):
        dir_from_file = np.loadtxt(datapath, delimiter=',', unpack=True)
        return [dir_mat.reshape((-1, model.dim)) for dir_mat in dir_from_file]

    def wrap_pca_dirs(model, dir_mat_list):
        return [GeneratedPCADirs(model, -1, -1, dir_mat=mat) for mat in dir_mat_list]

    def wrap_lin_dirs(model, dir_mat_list):
        return [GeneratedLinDirs(model, -1, -1, dir_mat=mat) for mat in dir_mat_list]


class Experiment(ABC):

    def __init__(self, *inputs, label=""):
        self.inputs = inputs
        self.plot = Plot()

        'Assuming that all models use the same dynamics and same initial set for now'
        self.model = inputs[0]['model']
        self.label = label

    """
    Execute experiment and dump results into spreadsheet.
    """
    @abstractmethod
    def execute(self, num_trials):
        pass

    """
    Assign directions based on pre-generated directions for each trial.
    """
    def assign_dirs(self, strat, trial_num, gen_dirs):
        if gen_dirs is not None:
            if isinstance(strat, MultiStrategy):
                for st in strat.strats:
                    self.__assign_dirs_by_strat(st, trial_num, gen_dirs)
            else:
                self.__assign_dirs_by_strat(strat, trial_num, gen_dirs)

    """
    Auxiliary method to assign directions based on strategy type.
    """
    def __assign_dirs_by_strat(self, strat, trial_num, gen_dirs):
        if isinstance(strat, AbstractPCAStrat):
            strat.pca_dirs = gen_dirs[trial_num][0]
        elif isinstance(strat, AbstractLinStrat):
            strat.lin_dirs = gen_dirs[trial_num][1]
        else:
            raise RuntimeError("Strategies have to be of either PCA, LinApp type.")

    """
    Method to load pre-generated directions from data directory. If not, pre-generate with supplied parameters and save to the data directory.
    model
    """
    def load_dirs(self, num_steps, num_trajs, num_trials):
        if KaaSettings.UsePreGenDirs:
            try:
                gen_dirs = DirSaveLoader.load_dirs(self.model, num_steps, num_trajs, KaaSettings.RandSeed)
                Output.prominent(f"Loaded directions from {KaaSettings.DataDir}")
            except IOError:
                Output.warning("WARNING: PRE-GENERATED DIRECTIONS NOT FOUND ON DISK. GENERATING DIRECTIONS FOR EXPERIMENT.")
                gen_dirs = self.__generate_dirs(num_steps, num_trajs, num_trials)
                Output.prominent("SAVING TO DISK.")
                DirSaveLoader.save_dirs(model, num_steps, num_trajs, KaaSettings.RandSeed, gen_dirs)
        else:
            gen_dirs = None

        return gen_dirs

    """
    Generate directions for each trial by incrementing random seed and generating both PCA and LinApp directions.
    """
    def __generate_dirs(self, num_steps, num_trajs, num_trials):
        generated_dirs = []
        for trial_num in range(num_trials):
            Output.prominent(f"GENERATED DIRECTIONS FOR TRIAL {trial_num} WITH {num_trajs} TRAJS FOR {num_steps} STEPS")
            ca_dirs = GeneratedPCADirs(self.model, num_steps, num_trajs)
            in_dirs = GeneratedLinDirs(self.model, num_steps, num_trajs)
            generated_dirs.append((gen_pca_dirs, gen_lin_dirs))
            update_seed()

        reset_seed()
        return generated_dirs

    """
    Saves data into a desired cell in spreadsheet.
    """
    def save_data_into_sheet(self, spreadsheet, trial_num, num_trials, flow_label, data):
        workbook = spreadsheet.workbook
        row_dict = spreadsheet.row_dict

        column_offset = trial_num
        row_offset = row_dict[flow_label]

        sheet = workbook.active
        sheet[chr(66 + column_offset) + str(row_offset)] = data

        if column_offset == num_trials - 1:
            sheet[chr(66 + num_trials) + str(row_offset)] = f"=AVERAGE(B{row_offset}:{chr(66 + num_trials - 1)}{row_offset})"
            sheet[chr(66 + num_trials + 1) + str(row_offset)] = f"=STDEV(B{row_offset}:{chr(66 + num_trials - 1)}{row_offset})"

        workbook.save(filename=os.path.join(PlotSettings.default_fig_path, self.label + '.xlsx'))

    """
    Initializes openpyxl spreadsheet to dump resulting data.
    """
    def generate_sheet(self, num_trials, row_labels=None):
        workbook = Workbook()
        sheet = workbook.active
        sheet.append(["Strategy"] + [f"Trial {i+1}" for i in range(num_trials)] + ["Mean", "Stdev"])

        'Initialize label-row dictionary'
        row_dict = {experi_input['label'] : row_idx + 2 for row_idx, experi_input in enumerate(self.inputs)}

        for experi_input in self.inputs:
            flow_label = experi_input['label']
            row = row_dict[flow_label]
            sheet['A' + str(row)] = flow_label

        workbook.save(filename=os.path.join(PlotSettings.default_fig_path, self.label + '.xlsx'))
        return SpreadSheet(workbook, row_dict)

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


    def gather_vol_data(self, experi_input):
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

"""
Experiment to compute the reachable set and estimate the total volume of all of its overapproximations.
"""
class VolumeExperiment(Experiment):

    def __init__(self, *inputs, label="Experiment"):
        super().__init__(*inputs, label=label)

    def execute(self, num_trials):
        num_steps = max([input['num_steps'] for input in self.inputs])
        num_trajs = max([input['num_trajs'] for input in self.inputs])

        spreadsheet = self.generate_sheet(num_trials)
        loaded_dirs = self.load_dirs(num_steps, num_trajs, num_trials)

        for experi_input in self.inputs:
            experi_strat = experi_input['strat']

            for trial_num in range(num_trials):
                Output.prominent(f"\n RUNNING EXPERIMENT {experi_input['label']} TRIAL:{trial_num} \n")
                self.assign_dirs(experi_strat, trial_num, loaded_dirs)

                flow_label, flow_vol = self.gather_vol_data(experi_input)
                self.save_data_into_sheet(spreadsheet, trial_num, num_trials, flow_label, flow_vol)

                experi_strat.reset()

"""
Experiment to measure deviations between generated directions for a strategy type over the course of the reachable set computation.
"""
class DeviationExperiment(Experiment):

    def __init__(self, *inputs, experi_type, label="Experiment"):
        super().__init__(*inputs, label=label)
        self.experi_type = experi_type

    def execute(self, num_trials):
        idx = dict(PCADev=0,
                   LinDev=1)

        strat_dirs_by_input = []
        row_labels = []
        for experi_input in self.inputs:
            num_steps = experi_input['num_steps']
            num_trajs = experi_input['num_trajs']
            label = experi_input['label']

            for dir_tuple in loaded_dirs:
                strat_dirs = dir_tuple[idx[self.experi_type]] #experi_type determines index of tuple to fetch (first index for PCA, second for LinApp)
                strat_dirs_by_input.append(pca_dirs)
                row_labels.append(label)

        for trial_num in range(num_trials):
            for row_label, strat_dirs_prev, strat_dirs_curr in zip(row_labels, strat_dirs_by_traj, strat_dirs_by_traj[1:]):
                prev_dirs = strat_dirs_prev[trial_num] #Corresponding directions from previous input to compare against
                curr_dirs = strat_dirs_curr[trial_num] #Corresponding directions from current input to compare against

                dirs_dist = self.__calc_dirs_dist(prev_dirs, curr_dirs)
                self.save_data_into_sheet(spreadsheet, trial_num, num_trials, row_label, dirs_dist)

    def __calc_dirs_dist(self, gen_dirs_one, gen_dirs_two):
         norm_dir_one = (gen_dirs_one.dir_mat.T / np.linalg.norm(gen_dirs_one.dir_mat, axis=1)).T
         norm_dir_two = (gen_dirs_two.dir_mat.T / np.linalg.norm(gen_dirs_two.dir_mat, axis=1)).T
         abs_dot_prods = np.abs(np.einsum('ij,ij->i', norm_dir_one, norm_dir_two))
         return np.min(abs_dot_prods)

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

"""
Update global random seed.
"""
def update_seed(offset=1):
    KaaSettings.RandSeed += offset

"""
Reset global random seed.
"""
def reset_seed():
    KaaSettings.RandSeed = 897987178

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
