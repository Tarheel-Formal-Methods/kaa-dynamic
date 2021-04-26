from collections import namedtuple
from itertools import product
from abc import ABC, abstractmethod
from pathlib import Path
import numpy as np
import os

from kaa.reach import ReachSet
from kaa.plotutil import *
from kaa.trajectory import Traj, TrajCollection
from kaa.settings import PlotSettings, KaaSettings
from kaa.templates import MultiStrategy, GeneratedDirs
from kaa.temp.pca_strat import AbstractPCAStrat, GeneratedPCADirs
from kaa.temp.lin_app_strat import AbstractLinStrat, GeneratedLinDirs
from kaa.log import Output
from kaa.timer import Timer

GenDirsTuple = namedtuple('GenDirsTuple', ['GenPCADirs', 'GenLinDirs'])

"""
Class responsible for saving/oading pre-generated directions from disk.
"""
class DirSaveLoader:

    """
    Loads pregenerated directions and sampled points used in calculating the directions from path specified
    in KaaSettings.DataDir. Returns a list of GenDirsTuples with initialized GeneratedDirs objects.
    @params model: Model object
            num_steps: Number of steps computed
            num_trajs: Number of trajectories used to compute directions
            seed: Random seed used to sample trajectory starting points in the initial box.

    @returns List of GenDirsTuples ordered by trial number
    """
    @staticmethod
    def load_dirs(model, num_steps, num_trajs, seed):
        pca_dirs_from_file = np.load(os.path.join(KaaSettings.DataDir, f"PCA{model}(T:{num_trajs})(Steps:{num_steps})(Seed:{seed}).npy"))
        lin_dirs_from_file = np.load(os.path.join(KaaSettings.DataDir, f"Lin{model}(T:{num_trajs})(Steps:{num_steps})(Seed:{seed}).npy"))
        samp_pts_from_file = np.load(os.path.join(KaaSettings.DataDir, f"SamPts{model}(T:{num_trajs})(Steps:{num_steps})(Seed:{seed}).npy"))

        pca_gendir_obj_list = DirSaveLoader.wrap_pca_dirs(model, pca_dirs_from_file, samp_pts_from_file)
        lin_gendir_obj_list = DirSaveLoader.wrap_lin_dirs(model, lin_dirs_from_file, samp_pts_from_file)

        gen_dirs_list = []
        for pca_gendir_obj, lin_gendir_obj in zip(pca_gendir_obj_list, lin_gendir_obj_list):
            gen_dirs_tuple = GenDirsTuple(pca_gendir_obj, lin_gendir_obj)
            gen_dirs_list.append(gen_dirs_tuple)

        return gen_dirs_list

    """
    Saves pregenerated directions to path specified by KaaSettings.DataDir. The files will be .npy files and their names will be designated according to
    the number of steps used in the direction computation, the number of trajectories used to compute the directions, and the random seed used to
    pick starting points within the initial box.
    @params model: Model object
            num_steps: Number of steps
            num_trajs: Number of trajectories
            seed: random seed
            gen_dirs_list: List of GenDirsTuple ordered by trial number (See above for definition.)
    """
    @staticmethod
    def save_dirs(model, num_steps, num_trajs, seed, gen_dirs_list):

        pca_dirs_by_trial = []
        lin_dirs_by_trial = []
        sampled_pts_by_trial = []

        for gen_dirs_tup in gen_dirs_list:
            gen_pca_dirs = gen_dirs_tup.GenPCADirs
            gen_lin_dirs = gen_dirs_tup.GenLinDirs

            pca_dirs_mat = gen_pca_dirs.dir_mat
            lin_dirs_mat = gen_lin_dirs.dir_mat

            sampled_points_mat = gen_pca_dirs.sampled_points

            pca_dirs_by_trial.append(pca_dirs_mat)
            lin_dirs_by_trial.append(lin_dirs_mat)
            sampled_pts_by_trial.append(sampled_points_mat)

        #print(np.asarray(sampled_pts_by_trial).shape)

        np.save(os.path.join(KaaSettings.DataDir, f"PCA{model}(T:{num_trajs})(Steps:{num_steps})(Seed:{seed}).npy"), pca_dirs_by_trial)
        np.save(os.path.join(KaaSettings.DataDir, f"Lin{model}(T:{num_trajs})(Steps:{num_steps})(Seed:{seed}).npy"), lin_dirs_by_trial)
        np.save(os.path.join(KaaSettings.DataDir, f"SamPts{model}(T:{num_trajs})(Steps:{num_steps})(Seed:{seed}).npy"), sampled_pts_by_trial)


    """
    Takes list of pregenerated PCA directions fetched from disk and initializes a GeneratedPCADirs
    object with the required starting data.
    @params model: Model object
            pca_dir_mat_list: List of PCA directions matrices ordered by trial. Each trial should have different random seed.
            sampled_pts_list: List of sampled point matrices.

    @returns List of GeneratedPCADirs objects ordered by trial number.
    """
    def wrap_pca_dirs(model, pca_dir_mat_list, sampled_pts_list):
        gen_pca_dirs_list = []

        #print(f"Sample Points Matrix List dim: sampled_pts_list.shape}")

        for mat_idx, pca_mat in enumerate(pca_dir_mat_list):
            sampled_pts_mat = sampled_pts_list[mat_idx]
            #print(f"Sample Points Matrix dim: {sampled_pts_mat.shape}")
            gen_pca_dirs = GeneratedPCADirs(model, -1, -1,
                                            dir_mat = pca_mat,
                                            sampled_points = sampled_pts_mat)

            gen_pca_dirs_list.append(gen_pca_dirs)

        return gen_pca_dirs_list


    """
    Takes list of pregenerated LinApp directions fetched from disk and initializes a GeneratedLinDirs
    object with the required starting data.
    @params model: Model object
            pca_lin_mat_list: List of LinApp directions matrices ordered by trial. Each trial should have different random seed.
            sampled_pts_list: List of sampled points

    @returns List of GeneratedLinDirs objects ordered by trial number.
    """
    def wrap_lin_dirs(model, lin_dir_mat_list, sampled_pts_list):
        gen_lin_dirs_list = []

        for mat_idx, lin_mat in enumerate(lin_dir_mat_list):
            sampled_pts_mat = sampled_pts_list[mat_idx]
            gen_lin_dirs = GeneratedLinDirs(model, -1, -1,
                                            dir_mat = lin_mat,
                                            sampled_points = sampled_pts_mat)

            gen_lin_dirs_list.append(gen_lin_dirs)

        return gen_lin_dirs_list

"""
Class containing all routines to conduct experiments with various strategies, extract
their output data, and plot/compare them.
"""
class Experiment(ABC):

    def __init__(self, *inputs, reach_comp_mode, label=""):
        self.inputs = inputs
        self.plot = Plot()

        'Assuming that all models use the same dynamics and same initial set for now'
        self.model = inputs[0]['model']
        self.label = label
        self.reach_comp_mode = reach_comp_mode

    """
    Execute experiment and dump results into spreadsheet.
    """
    @abstractmethod
    def execute(self, num_trials):
        pass

    @property
    def max_num_steps(self):
        return max([experi_input['num_steps'] for experi_input in self.inputs])

    @property
    def max_num_trajs(self):
        return max([experi_input['num_trajs'] for experi_input in self.inputs])

    """
    Assign directions based on pre-generated directions for each trial.
    @params strat: TempStrategy object
            trial_num: Trial number
            gen_dirs_list: List of GenDirsTuples ordered by trial number.
    """
    def assign_dirs(self, strat, trial_num, gen_dirs_list):
        if gen_dirs_list is not None:
            if isinstance(strat, MultiStrategy):
                for st in strat:
                    self.__assign_dirs_by_strat(st, trial_num, gen_dirs_list)
            else:
                self.__assign_dirs_by_strat(strat, trial_num, gen_dirs_list)

    """
    Auxiliary method to assign directions based on strategy type.
    @params strat: TempStrategy object
            trial_num: Trial number
            gen_dirs_list: List of GenDirsTuples ordered by trial number
    """
    def __assign_dirs_by_strat(self, strat, trial_num, gen_dirs_list):
        if isinstance(strat, AbstractPCAStrat):
            strat.dirs = gen_dirs_list[trial_num].GenPCADirs

        elif isinstance(strat, AbstractLinStrat):
            strat.dirs = gen_dirs_list[trial_num].GenLinDirs

        else:
            raise RuntimeError("Strategies have to be of either PCA, LinApp type.")

    """
    Method to load pre-generated directions from data directory. If not, pre-generate with supplied parameters and save to the data directory.
    @params num_steps: Number of steps to be used for computation.
            num_trajs: Number of trajectories used during computation.
            seed: Random seed used for computation.

    @returns List of GenDirsTuples ordered by trial number.
    """
    def load_dirs(self, num_steps, num_trajs, num_trials):
        try:
            gen_dirs_list = DirSaveLoader.load_dirs(self.model, num_steps, num_trajs, KaaSettings.RandSeed)
            Output.prominent(f"Loaded directions from {KaaSettings.DataDir}")

        except IOError:
            Output.warning("WARNING: PRE-GENERATED DIRECTIONS NOT FOUND ON DISK. GENERATING DIRECTIONS FOR EXPERIMENT.")
            gen_dirs_list = self.__generate_dirs(num_steps, num_trajs, num_trials)

            Output.prominent("SAVING TO DISK.")
            DirSaveLoader.save_dirs(self.model, num_steps, num_trajs, KaaSettings.RandSeed, gen_dirs_list)

        return gen_dirs_list

    """
    Generate directions for each trial by incrementing random seed and generating both PCA and LinApp directions.
    @params num_steps: Number of steps to propagate trajectories according to system dynamics
            num_trajs: Number of trajectories to use to generate directions.
            num_trials: Number of trials to increment seed and generate new set of directions.

    @returns List of GenDirsTuples ordered by trial number.
    """
    def __generate_dirs(self, num_steps, num_trajs, num_trials):
        generated_dirs = []

        for trial_num in range(num_trials):
            Output.prominent(f"GENERATED DIRECTIONS FOR TRIAL {trial_num} WITH {num_trajs} TRAJS FOR {num_steps} STEPS")
            gen_pca_dirs = GeneratedPCADirs(self.model, num_steps, num_trajs)
            gen_lin_dirs = GeneratedLinDirs(self.model, num_steps, num_trajs)

            gen_dirs_tuple = GenDirsTuple(gen_pca_dirs, gen_lin_dirs)
            generated_dirs.append(gen_dirs_tuple)

            update_seed()

        reset_seed()
        return generated_dirs

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

    def calc_flowpipe(self, experi_input):
        try:
            model = experi_input['model']
            strat = experi_input['strat']
            flow_label = experi_input['label']
            num_steps = experi_input['num_steps']

            mod_reach = ReachSet(model, strat=strat, label=flow_label, reachcompmode=self.reach_comp_mode)
            mod_flow = mod_reach.computeReachSet(num_steps)

            return mod_flow

        except Exception as excep:
            raise
            return (experi_input['label'], str(excep))

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
    Sample points from the edges of the box and propagate them for a number of steps.
    """
    def simulate_border_points(self, num_steps):
        init_box_inter = self.model.init_box
        border_points = self.__get_init_box_borders(init_box_inter)

        trajs = [Traj(self.model, point, num_steps) for point in border_points]

        return TrajCollection(self.model, trajs)

    """
    Find corner vertices for an initial box along with midpoints between the corners.
    @params init_box : intervals of the box given as a list of lists where each member's left,right value
                       are the start,end points respectively for the intervals of the box.
    @returns list of border points.
    """
    def __get_init_box_borders(self, init_box):
        midpoints = [start + (end - start) / 2 for start, end in init_box]
        border_points = list(product(*init_box))

        for point_idx, point in enumerate(midpoints):
            half_points = [init_inter if point_idx != inter_idx else [point] for inter_idx, init_inter in enumerate(init_box)]
            border_points += list(product(*half_points))

        return border_points

    """
    Initialize strategies based on input paramters supplied to Experiment object. Flags and some attributes
    must be initialized for the strategies to function properly during the reachable set computation.
    @params experi_input: dictionary containing experi_input
            num_trials: number of trials for experiment
    """
    def initialize_strat(self, experi_input, num_trials):
        experi_strat = experi_input['strat']
        loaded_dirs = None

        if not experi_strat: return #If strat is None, nothing to initialize

        experi_supp_mode = experi_input['supp_mode']
        experi_pregen_mode = experi_input['pregen_mode']

        if experi_supp_mode:
            if isinstance(experi_strat, MultiStrategy):
                for strat in experi_strat:
                    strat.use_supp_points = True
            else:
                experi_strat.use_supp_points = True

        elif experi_pregen_mode:
            experi_num_trajs = experi_input['num_trajs']
            experi_num_steps = experi_input['max_steps']

            experi_strat.num_trajs =  experi_num_trajs

            loaded_dirs = self.load_dirs(experi_num_steps, experi_num_trajs, num_trials)
            self.assign_dirs(experi_strat, 0, loaded_dirs)

        return loaded_dirs

    def print_input_params(self, experi_input, trial_num=-1):
        experi_strat = experi_input['strat']
        Output.prominent(f"Running Experiment {experi_input['label']} Trial:{trial_num if trial_num >= 0 else None}")

        if experi_strat:
            experi_supp_mode = experi_input['supp_mode']
            experi_pregen_mode = experi_input['pregen_mode']
            experi_num_trajs = experi_input['num_trajs']

            Output.prominent(f"Using following parameters for experiments:")
            Output.prominent(f"Use Support Points: {experi_supp_mode}")
            Output.prominent(f"Use Pre-gen Points: {experi_pregen_mode}")

            if experi_pregen_mode:
                Output.prominent(f"Number of Trajectories Used: {experi_num_trajs}")

    def __str__(self):
        return self.label

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
