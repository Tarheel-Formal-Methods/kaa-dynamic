import numpy as np

from kaa.trajectory import Traj, TrajCollection
from kaa.log import Output
from kaa.experiment import Experiment
from kaa.bundle import BundleTransMode
from kaa.temp.random_static_strat import RandomStaticStrat
from kaa.temp.diag_static_strat import RandomDiagStaticStrat
from kaa.temp.lin_app_strat import SlidingLinStrat
from kaa.temp.pca_strat import SlidingPCAStrat
from kaa.templates import MultiStrategy
from kaa.spreadsheet import *
from kaa.timer import Timer
from kaa.dataload import CovidDataLoader
from models.covid import Covid_UnitBox, Covid
from kaa.linearsystem import calc_box_volume

"""
Experiment to generate statistics useful for profiling Kaa
"""


class ProfileExperiment(Experiment):
    def __init__(self, *inputs, label="Experiment"):
        super().__init__(*inputs, label=label)

    def execute(self):
        for experi_input in self.inputs:
            self.print_input_params(experi_input)

            self.initialize_strat(experi_input, 10)
            self.calc_flowpipe(experi_input)

        Timer.generate_stats()


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
                strat_dirs = dir_tuple[idx[
                    self.experi_type]]  # experi_type determines index of tuple to fetch (first index for PCA, second for LinApp)
                strat_dirs_by_input.append(pca_dirs)
                row_labels.append(label)

        spreadsheet = SpreadSheet(self.model, self.label)
        spreadsheet.generate_sheet(self.inputs, num_trials)

        for trial_num in range(num_trials):
            for row_label, strat_dirs_prev, strat_dirs_curr in zip(row_labels, strat_dirs_by_traj,
                                                                   strat_dirs_by_traj[1:]):
                prev_dirs = strat_dirs_prev[
                    trial_num]  # Corresponding directions from previous input to compare against
                curr_dirs = strat_dirs_curr[trial_num]  # Corresponding directions from current input to compare against

                dirs_dist = self.__calc_dirs_dist(prev_dirs, curr_dirs)
                spreadsheet.save_data_into_sheet(trial_num, num_trials, row_label, dirs_dist)

    def __calc_dirs_dist(self, gen_dirs_one, gen_dirs_two):
        norm_dir_one = (gen_dirs_one.dir_mat.T / np.linalg.norm(gen_dirs_one.dir_mat, axis=1)).T
        norm_dir_two = (gen_dirs_two.dir_mat.T / np.linalg.norm(gen_dirs_two.dir_mat, axis=1)).T
        abs_dot_prods = np.abs(np.einsum('ij,ij->i', norm_dir_one, norm_dir_two))
        return np.min(abs_dot_prods)


"""
Experiment to calculate and plot the phase plot.
"""


class PhasePlotExperiment(Experiment):

    def __init__(self, *inputs, separate=False, plot_border_traj=False):
        super().__init__(*inputs)
        self.separate = separate
        self.plot_border_traj = plot_border_traj

    def populate_plot(self):
        num_steps = self.max_num_steps

        if self.plot_border_traj:
            self.plot.add(self.simulate_border_points(num_steps))

        for experi_input in self.inputs:
            self.initialize_strat(experi_input, 10)
            self.plot.add(self.calc_flowpipe(experi_input))

    def plot_flowpipes(self, *var_tup, xlims, ylims):
        self.plot.plot({'type': 'Phase',
                        'vars': var_tup,
                        'separate_flag': self.separate,
                        'xlims': xlims,
                        'ylims': ylims})

    def execute(self, *var_tup, xlims=None, ylims=None):
        self.populate_plot()
        self.plot_flowpipes(*var_tup, xlims=xlims, ylims=ylims)


class InitReachPlotExperiment(Experiment):

    def __init__(self, *inputs, log_scale=False):
        super().__init__(*inputs)
        self.log_scale_flag = log_scale

    def execute(self):
        for experi_input in self.inputs:
            self.initialize_strat(experi_input, 10)  # fix this
            # self.print_input_params(experi_input, trial_num=None)
            self.plot.add(self.calc_flowpipe(experi_input))

        self.plot.plot({'type': 'InitVolReachVol',
                        'flowpipe_indepen_data': None,
                        'log_scale_flag': self.log_scale_flag})


class InitReachVSRandomPlotExperiment(Experiment):
    def __init__(self, *inputs, num_ran_temps=20, num_trials=10, ran_diag=False, log_scale=False, precalc_vals=None):
        super().__init__(*inputs)
        self.num_trials = num_trials
        self.num_ran_temps = num_ran_temps
        self.ran_diag_flag = ran_diag
        self.log_scale_flag = log_scale
        self.precalc_vals = precalc_vals

    def execute(self):
        ran_flow_data_tups = []
        for input_idx, experi_input in enumerate(self.inputs):
            self.initialize_strat(experi_input, 10)
            self.print_input_params(experi_input)

            ran_data_tup = self.__avg_ran_flowpipes(experi_input, input_idx)
            ran_flow_data_tups.append(ran_data_tup)
            self.plot.add(self.calc_flowpipe(experi_input))

        self.plot.plot({'type': 'InitVolReachVol',
                        'flowpipe_indepen_data': ran_flow_data_tups,
                        'log_scale_flag': self.log_scale_flag})

    def __avg_ran_flowpipes(self, experi_input, input_idx):
        input_model = experi_input['model']
        input_traj = experi_input['num_trajs']
        input_supp = experi_input['supp_mode']
        input_pregen = experi_input['pregen_mode']
        num_steps = experi_input['num_steps']

        ran_label = f"{input_model.name} Random Static {'Diag' if self.ran_diag_flag else ''} Strat 10 trials"

        if self.precalc_vals:
            assert len(self.precalc_vals) == len(
                self.inputs), "Averaged precalculated flowpipe values must correspond to each input."
            return (ran_label, calc_box_volume(input_model.init_box), self.precalc_vals[input_idx])

        trial_vol_arr = np.empty(self.num_trials)

        for trial in range(self.num_trials):
            ran_strat = (RandomDiagStaticStrat(input_model, self.num_ran_temps) if self.ran_diag_flag else
                         RandomStaticStrat(input_model, self.num_ran_temps))

            ran_experi_input = dict(model=input_model,
                                    strat=ran_strat,
                                    label=ran_label,
                                    supp_mode=input_supp,
                                    pregen_mode=input_pregen,
                                    num_trajs=input_traj,
                                    num_steps=num_steps)

            self.print_input_params(ran_experi_input, trial_num=trial)
            trial_vol_arr[trial] = self.calc_flowpipe(ran_experi_input).total_volume
            # spreadsheet.save_data_into_sheet(ran_label, trial_num, (flow_vol, f"{reach_time[0]} min {reach_time[1]} sec)"))

            Output.write(f"INPUT: {input_model.init_box} Trial {trial}: {trial_vol_arr[trial]}")

        return (ran_label, calc_box_volume(input_model.init_box), np.average(trial_vol_arr))


class VolumeDataExperiment(Experiment):
    def __init__(self, *inputs, label="VolumeDataExperiment", num_trials=1):
        super().__init__(*inputs, label=label)
        self.num_trials = num_trials

    def execute(self):
        num_steps = self.max_num_steps

        spreadsheet = SpreadSheet(self.model, self.label)
        spreadsheet.generate_sheet(self.inputs, ("Volume", "Time"), self.num_trials)

        for experi_input in self.inputs:
            self.print_input_params(experi_input)
            self.initialize_strat(experi_input, 10)
            experi_strat = experi_input['strat']

            for trial_num in range(self.num_trials):
                flowpipe = self.calc_flowpipe(experi_input)

                flow_label = experi_input['label']
                flow_vol = flowpipe.total_volume
                reach_time = divmod(Timer.generate_stats(), 60)

                spreadsheet.save_data_into_sheet(flow_label, trial_num,
                                                 (flow_vol, f"{reach_time[0]} min {reach_time[1]} sec)"))
                experi_strat.reset()


class VolumePlotExperiment(Experiment):

    def __init__(self, *inputs, label="VolumePlotExperiment", accum=True, plot_all_vol=False):
        super().__init__(*inputs, label=label)
        self.accum_flag = accum
        self.plot_all_vol_flag = plot_all_vol

    def execute(self):
        num_steps = self.max_num_steps

        spreadsheet = SpreadSheet(self.model, self.label)
        spreadsheet.generate_sheet(self.inputs, ("Volume", "Time"), 1)

        for experi_input in self.inputs:
            self.print_input_params(experi_input)
            self.initialize_strat(experi_input, 10)

            flowpipe = self.calc_flowpipe(experi_input)

            flow_label = experi_input['label']
            flow_vol = flowpipe.total_volume
            reach_time = divmod(Timer.generate_stats(), 60)

            spreadsheet.save_data_into_sheet(flow_label, 0, (flow_vol, f"{reach_time[0]} min {reach_time[1]} sec)"))
            self.plot.add(flowpipe)

        self.plot.plot({'type': 'Volume',
                        'accum_flag': self.accum_flag,
                        'plot_all_vol_flag': self.plot_all_vol_flag})


"""
Experiment to calculate and plot the projection reachable sets.
"""


class ProjectionPlotExperiment(Experiment):

    def __init__(self, *inputs, separate=False, plot_border_traj=False, plot_total_width=False):
        super().__init__(*inputs)
        self.separate_flag = separate
        self.plot_border_traj = plot_border_traj
        self.plot_total_width_flag = plot_total_width

    def populate_plot(self):
        num_steps = self.max_num_steps

        if self.plot_border_traj:
            self.plot.add(self.simulate_border_points(num_steps))

        for experi_input in self.inputs:
            self.print_input_params(experi_input)

            self.initialize_strat(experi_input, 10)
            self.plot.add(self.calc_flowpipe(experi_input))

    def plot_flowpipes(self, *var_tup, xlims, ylims):
        self.plot.plot({'type': 'Projection',
                        'vars': var_tup,
                        'plot_width_flag': self.plot_total_width_flag,
                        'xlims': xlims,
                        'ylims': ylims})

    def execute(self, *var_tup, xlims=None, ylims=None):
        self.populate_plot()
        self.plot_flowpipes(*var_tup, xlims=xlims, ylims=ylims)


"""
Experiment created to plot Covid models against real-world data for Indian Supermodel
(https://www.iith.ac.in/~m_vidyasagar/arXiv/Super-Model.pdf)
"""


class CovidDataPlotExperiment(Experiment):

    def __init__(self, earliest_date=None, latest_date=None, filename='case_time_series.csv', total_pop=1.36E9,
                 epsilon=1E-6):
        self.filename = filename

        data_loader = CovidDataLoader(self.filename)
        self.data_dict = data_loader.fetch_data(earliest_date, latest_date)

        self.earliest_date = earliest_date if earliest_date else self.data_dict.items()[0][0]
        self.latest_date = latest_date if earliest_date else self.data_dict.items()[-1][0]
        self.total_pop = total_pop
        self.epsilon = epsilon

        input = self.__init_inputs()
        super().__init__(input)
        self.__calc_trajs()
        self.input = self.inputs[0]

    def execute(self):
        self.print_input_params(self.input)
        self.initialize_strat(self.input, 10)

        self.plot.add(self.calc_flowpipe(self.input))

        self.plot.plot({'type': 'CovidProj',
                        'data_dict': self.data_dict,
                        'steps_in_day': self.num_steps_in_day,
                        'total_pop': self.total_pop})

    def __init_inputs(self):
        init_params = self.__estimate_init_params(self.total_pop)
        num_days = len(self.data_dict) - 1

        confirmed = init_params['Confirmed']
        recovered = init_params['Recovered']
        deceased = init_params['Deceased']

        '''
        According to Sutra paper (https://arxiv.org/abs/2101.09158), 
        asymptomatic patients tend to be 80% of all '
        tested positive. 
        '''
        self.init_confirmed_asymp = 0.8 * confirmed
        self.init_confirmed_symp = 0.2 * confirmed

        self.init_recovered_asymp = 0.8 * recovered
        self.init_recovered_symp = 0.2 * recovered

        init_suspectible = 1 - confirmed - recovered - deceased

        self.init_suspectible_asymp = 0.8 * init_suspectible
        self.init_suspectible_symp = 0.2 * init_suspectible

        self.init_deceased = deceased

        init_box = ((self.init_suspectible_asymp - self.epsilon, self.init_suspectible_asymp + self.epsilon),
                    (self.init_suspectible_symp - self.epsilon, self.init_suspectible_symp + self.epsilon),
                    (self.init_confirmed_asymp - self.epsilon, self.init_confirmed_asymp + self.epsilon),
                    (self.init_confirmed_symp - self.epsilon, self.init_confirmed_symp + self.epsilon),
                    (self.init_recovered_asymp - self.epsilon, self.init_recovered_asymp + self.epsilon),
                    (self.init_recovered_symp - self.epsilon, self.init_recovered_symp + self.epsilon),
                    (self.init_deceased - self.epsilon, self.init_deceased + self.epsilon),
                    (0.108, 0.112),
                    (0.078, 0.082)
                    )

        self.model = Covid_UnitBox(delta=0.5, init_box=init_box)

        self.num_steps_in_day = int(1 / self.model.step_size)
        self.total_num_steps = self.num_steps_in_day * num_days - 1

        pca_window_size = 8
        lin_window_size = 0

        pca_strat = SlidingPCAStrat(self.model, lifespan=pca_window_size)
        lin_strat = SlidingLinStrat(self.model, lifespan=lin_window_size)

        input = dict(model=self.model,
                     strat=MultiStrategy(pca_strat, lin_strat),
                     label=f"Covid PCA WinSize {pca_window_size} and Lin WinSize {lin_window_size}",
                     supp_mode=True,
                     pregen_mode=False,
                     num_trajs=0,
                     num_steps=self.total_num_steps,
                     trans_mode=BundleTransMode.AFO)

        return input

    def __calc_trajs(self):
        beta_step_size = 0.002
        traj_list = []
        for beta in np.arange(0.108, 0.117, beta_step_size):
            beta_traj = Traj(self.model,
                             (self.init_suspectible_asymp,
                              self.init_suspectible_symp,
                              self.init_confirmed_asymp,
                              self.init_confirmed_symp,
                              self.init_recovered_asymp,
                              self.init_recovered_symp,
                              self.init_deceased,
                              beta,
                              0.08
                              ),
                             self.total_num_steps,
                             label=f"β = {beta}, γ = 0.08")
            traj_list.append(beta_traj)

        self.plot.add(TrajCollection(self.model, traj_list))

    def __estimate_init_params(self, total_pop):
        initial_figures = self.data_dict[self.earliest_date]

        total_confirmed = int(initial_figures['Confirmed'])
        total_recovered = int(initial_figures['Recovered'])
        total_death = int(initial_figures['Deceased'])

        init_params_dict = {'Confirmed': total_confirmed / total_pop,
                            'Recovered': total_recovered / total_pop,
                            'Deceased': total_death / total_pop}

        return init_params_dict


class CompAniExperiment(Experiment):

    def __init__(self, *inputs):
        super().__init__(*inputs)

    def execute(self, x, y, ptope_order, filename, plot_pts=None):
        if not plot_pts: plot_pts = [False for _ in enumerate(self.inputs)]

        flowpipes = []
        for experi_input in self.inputs:
            self.initialize_strat(experi_input, 10)
            flowpipes.append(self.calc_flowpipe(experi_input))

        animation = SlideCompareAnimation(*flowpipes)

        border_trajs = self.simulate_border_points(self.max_num_steps)
        animation.add(border_trajs)

        animation.animate(x, y, ptope_order, plot_pts, filename)
