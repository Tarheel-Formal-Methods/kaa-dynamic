import os
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as pat
import matplotlib.animation as animate
from matplotlib.ticker import MultipleLocator
from matplotlib.legend import Legend
from pathlib import Path
from datetime import date
from itertools import groupby
from abc import ABC, abstractmethod

from kaa.settings import PlotSettings
from kaa.trajectory import TrajCollection, Traj
from kaa.flowpipe import FlowPipe
from kaa.timer import Timer

plt.rcParams.update({'font.size': PlotSettings.PlotFont})


def plot_traj_proj(model, ax, trajs, x, num_steps, linewidth=1):
    if not trajs: return

    t = np.arange(0, num_steps, 1)

    for traj_idx, traj in enumerate(trajs):
        y_coord = traj.get_proj(x)
        y_coord = y_coord[:num_steps]

        if traj.label:
            ax.plot(t, y_coord, color=f"C{traj_idx + model.dim}", linewidth=linewidth, label=traj.label)
        else:
            ax.plot(t, y_coord, color=f"C{traj_idx + model.dim}", linewidth=linewidth)


def plot_trajs(model, ax, trajs, x, y, num_steps=0):
    if not trajs: return

    for traj_idx, traj in enumerate(trajs):
        x_coord = traj.get_proj(x)
        y_coord = traj.get_proj(y)

        if num_steps:
            x_coord = x_coord[:num_steps]
            y_coord = y_coord[:num_steps]

        ax.plot(x_coord, y_coord, color="k", linewidth=0.3)


class Subplot:

    def __init__(self, model, label, flowpipes, num_plots, num_steps):
        self.label = label
        self.num_plots = num_plots
        self.model = model
        self.num_steps = num_steps
        self.flowpipes = flowpipes
        self.num_flowpipes = len(flowpipes)

    @abstractmethod
    def plot(self):
        pass

    def print_final_vol_stats(self):
        for flowpipe in self.flowpipes:
            print(f"VOLUME OF FINAL SYSTEM FOR FLOWPIPE {flowpipe.label}: {flowpipe[-1].volume[1]}")


# Some sort of organization through inheritance is warranted.
class ProjectionSubplot(Subplot):

    def __init__(self, model, vars, flowpipes, num_steps, trajs, plot_width_flag, x_lims, y_lims, scale_factor):
        self.vars = vars
        self.trajs = trajs
        self.plot_width_flag = plot_width_flag
        self.x_lims = x_lims
        self.y_lims = y_lims
        self.scale_factor = scale_factor
        super().__init__(model,
                         "Projection",
                         flowpipes,
                         len(self.vars) + 1 if plot_width_flag else len(self.vars),
                         num_steps)

    """
    Plots the projections of the trajectories and flowpipes stored in Plot object.
    @params: *var_tup: indices of desired variables
              path: optional variable designated the path to store the generated matplotlib figure.
    """

    def plot(self, axs):
        num_vars = len(self.vars)
        num_flowpipes = len(self.flowpipes)
        name = self.model.name

        assert self.model is not None, "No data has been added to the Plot object."

        for ax_idx, (var_idx, ax) in enumerate(zip(self.vars, axs)):
            var = self.model.vars[var_idx]

            plot_traj_proj(self.model, ax, self.trajs, var_idx, self.num_steps)

            for flow_idx, flowpipe in enumerate(self.flowpipes):
                flow_min, flow_max = flowpipe.get_proj(var_idx)
                flowpipe_label = flowpipe.label

                t = np.arange(0, len(flowpipe), 1)

                if self.scale_factor:
                    t = t * (self.scale_factor / (len(flowpipe) - 1)) #For some reaosn, t*= gives type casting error.

                ax.fill_between(t, flow_min, flow_max,
                                label=flowpipe_label,
                                color=f"C{flow_idx}",
                                alpha=0.8 - (flow_idx / 10),
                                zorder=num_flowpipes - flow_idx)

            ax.set_xlabel("t: time steps", fontsize=PlotSettings.PlotFont)
            ax.set_ylabel(f"Reachable Set for {var}", fontsize=PlotSettings.PlotFont)
            ax.set_title(f"Projection of Reachable Set for {name} Variable: {var}")
            ax.legend(loc=2)

            ax.tick_params(axis='both', labelsize=PlotSettings.PlotFont)

            if self.x_lims:
                ax.set_xlim(self.x_lims)
            if self.y_lims:
                ax.set_ylim(self.y_lims)

        'Total width is displayed on the last subplot by default'
        if self.plot_width_flag:
            ax = axs[num_vars]

            for flow_idx, flowpipe in enumerate(self.flowpipes):
                flow_min, flow_max = flowpipe.get_total_width_reachable_set()
                flowpipe_label = flowpipe.label

                t = np.arange(0, len(flowpipe), 1)

                if self.scale_factor:
                    t *= self.scale_factor

                ax.fill_between(t, flow_min, flow_max,
                                label=flowpipe_label,
                                color=f"C{flow_idx}", alpha=0.4 + (flow_idx / 10),
                                zorder=num_flowpipes - flow_idx)

            ax.set_xlabel("t: time steps", fontsize=PlotSettings.PlotFont)
            ax.set_ylabel(f"Reachable Set for s", fontsize=PlotSettings.PlotFont)
            ax.set_title(f"Projection of Reachable Set for s")
            ax.legend(loc=2)

            ax.tick_params(axis='both', labelsize=PlotSettings.PlotFont)

            if self.x_lims:
                ax.set_xlim(self.x_lims)
            if self.y_lims:
                ax.set_ylim(self.y_lims)

"""
Assumes the dynamics shown in https://www.iith.ac.in/~m_vidyasagar/arXiv/Super-Model.pdf
"""


class CovidProjectionDateSubplot(Subplot):

    def __init__(self, model, flowpipes, num_steps, trajs, data_dict, steps_in_day, total_pop, time_interval_str,
                 category=None):
        self.data_dict = data_dict
        self.steps_in_day = steps_in_day
        self.total_pop = total_pop
        self.trajs = trajs
        self.time_interval_str = time_interval_str
        self.category = category
        super().__init__(model,
                         "CovidProj",
                         flowpipes,
                         1 if category else 2,
                         num_steps)

    def plot(self, axs):
        num_flowpipes = len(self.flowpipes)
        name = self.model.name
        step_size = self.model.step_size

        pos_dict = {'Confirmed': 1,
                    'Recovered': 2,
                    'Deceased': 3
                    }

        category_iter = [self.category] if self.category else ('Confirmed', 'Recovered')

        for ax, category in zip(axs, category_iter):
            assert category in set(pos_dict.keys()), f"{category} not a valid category in COVID dynamics"

            ax_idx = pos_dict[category]

            param_var_trajs = self.trajs.merge_dims(2, 3).scale(
                self.total_pop)  # Merge asymp and symp confirmed dimensions for trajectories
            plot_traj_proj(self.model,
                           ax,
                           param_var_trajs,
                           2,
                           self.num_steps,
                           linewidth=1.5)

            ax_data = [date_data[category] for date_data in self.data_dict.values()]
            num_of_days = self.num_steps * step_size

            for flow_idx, flowpipe in enumerate(self.flowpipes):
                flowpipe_label = flowpipe.label

                if self.category != 'Deceased':
                    first_flow_min, first_flow_max = flowpipe.get_proj(2 * pos_dict[category])
                    second_flow_min, second_flow_max = flowpipe.get_proj(2 * pos_dict[category] + 1)

                    flow_min = first_flow_min + second_flow_min
                    flow_max = first_flow_max + second_flow_max

                else:
                    flow_min, flow_max = flowpipe.get_proj(2 * pos_dict[category])

                scaled_flow_min = flow_min * self.total_pop
                scaled_flow_max = flow_max * self.total_pop

                t = np.arange(0, len(flowpipe), 1)
                ax.fill_between(t, scaled_flow_min, scaled_flow_max,
                                label="Predicted Region",
                                color=f"C{flow_idx}",
                                alpha=0.35 - (flow_idx / 10),
                                zorder=1)

            ax.plot(range(0, (self.num_steps - 1) + self.steps_in_day, self.steps_in_day),
                    ax_data,
                    'o-r',
                    label="Real Data",
                    linewidth=3.5)

            dates = list(self.data_dict.keys())
            day_tick_gap = 5
            date_labels = []
            for idx in range(self.num_steps + 1):
                wrapped_idx, remainder = divmod(idx, self.steps_in_day)
                date_labels.append(dates[wrapped_idx]
                                   if not remainder and (wrapped_idx % day_tick_gap == 0)
                                   else "")

            ax.set_xticks(range(self.num_steps + 1))
            ax.set_xticklabels(date_labels)
            ax.xaxis.set_major_locator(MultipleLocator(day_tick_gap))
            ax.xaxis.set_minor_locator(MultipleLocator(self.steps_in_day))
            ax.tick_params(which="major", width=6)

            ax.set_xlabel("Dates", fontsize=PlotSettings.PlotFont)
            ax.set_ylabel(f"Population for {category}", fontsize=PlotSettings.PlotFont)
            ax.set_title(f"Projection for {category} Against Real-World India Data for {self.time_interval_str}")
            ax.legend(loc=2)
            ax.tick_params(axis='both', labelsize=PlotSettings.PlotFont)


class PhaseSubplot(Subplot):

    def __init__(self, model, vars, flowpipes, num_steps, separate_flag, trajs, x_lims, y_lims):
        assert len(vars) == 2, "Phase is 2D, so vars can only be a 2-element tuple."
        self.x = vars[0]
        self.y = vars[1]
        self.separate_flag = separate_flag
        self.trajs = trajs
        self.x_lims = x_lims
        self.y_lims = y_lims
        super().__init__(model, "Phase", flowpipes, 1, num_steps)

    """
    Plots phase between two variables of dynamical system.
    NOTE: This method for now creates rather crude, segmented phase plots by simply scattering the support points.

    @params x: index of variable to be plotted as x-axis of desired phase
            y: index of variable to be plotted as y-axis of desired phase
    """

    def plot(self, ax):
        assert len(ax) == 1, "Only one axis object for plotting phase."

        Timer.start('Phase')
        ax = ax[0]
        for flow_idx, flowpipe in enumerate(self.flowpipes):
            self.__halfspace_inter_plot(ax, self.x, self.y, flowpipe, flow_idx, 0.8 - (flow_idx / 10))
            # self.__support_plot(flowpipe, flow_idx,  x, y, ax)

        plot_trajs(self.model, ax, self.trajs, self.x, self.y)
        self.__phase_plot_legend(self.x, self.y, ax)

        if self.x_lims:
            ax.set_xlim(self.x_lims)
        if self.y_lims:
            ax.set_ylim(self.y_lims)

        phase_time = Timer.stop('Phase')
        x_var, y_var = self.model.vars[self.x], self.model.vars[self.y]
        print("Plotting phase for dimensions {}, {} done -- Time Spent: {}".format(x_var, y_var, phase_time))

    """
    Use scipy.HalfspaceIntersection to fill in phase plot projections.
    @params: flowpipe: FlowPipe object to plot.
    """

    def __halfspace_inter_plot(self, ax, x, y, flowpipe, flow_idx, alpha):
        for bund in flowpipe:
            if not self.separate_flag:
                'Temp patch. Revise to start using Line, fontsize=25arSystems for future work.'
                self.__plot_halfspace(x, y, ax,
                                      bund.getIntersect(),
                                      alpha,
                                      flow_idx)
            else:
                for ptope_idx, ptope in enumerate(bund.all_ptopes):
                    self.__plot_halfspace(x, y, ax,
                                          ptope,
                                          alpha,
                                          flow_idx + ptope_idx)

    """
    Plot linear system through scipy.HalfspaceIntersection
    """

    def __plot_halfspace(self, x, y, ax, sys, alpha, idx_offset):
        'Routines in Bundle give None values for queries of non-existent ptopes'
        if sys is None: return sys

        dim = self.model.dim
        comple_dim = np.asarray([True if i in [x, y] else False for i in range(dim)])

        vertices = sys.vertices

        if vertices is None: return sys

        proj_vertices = np.unique(vertices[:, comple_dim], axis=0).tolist()
        center_pt = sys.chebyshev_center.center

        'Sort by polar coordinates to ensure proper plotting of boundary'
        proj_vertices.sort(key=lambda v: math.atan2(v[1] - center_pt[1], v[0] - center_pt[0]))

        ptope = pat.Polygon(proj_vertices, fill=True, color=f"C{idx_offset}", alpha=alpha)
        ax.add_patch(ptope)

        inter_x, inter_y = zip(*proj_vertices)
        ax.scatter(inter_x, inter_y, s=0.01)

    """
    Routine to populate the legend information for a phase plot.
    @params x: index of x variable
            y: index of y variable
            phase_ax: Axis object associated to phase plot.
            lims: optional data denoting x,y axis dimensions
    """

    def __phase_plot_legend(self, x, y, phase_ax):
        x_var, y_var = self.model.vars[x], self.model.vars[y]

        phase_ax.set_xlabel(f'{x_var}')
        phase_ax.set_ylabel(f'{y_var}')
        phase_ax.set_title("Projection of Phase Plot for {} Variables: {}".format(self.model.name, (x_var, y_var)))

        phase_ax.tick_params(axis='both', labelsize=PlotSettings.PlotFont)

        axis_patches = []
        for flow_idx, flowpipe in enumerate(self.flowpipes):
            axis_patches.append(pat.Patch(color=f"C{flow_idx}", label=flowpipe.label))

        phase_ax.legend(handles=axis_patches)

    def __support_plot(self, flowpipe, flow_idx, x, y, ax):
        dim = self.model.dim

        'Define the following projected normal vectors.'
        norm_vecs = np.zeros([4, dim])
        norm_vecs[0][x] = 1;
        norm_vecs[1][y] = 1;
        norm_vecs[2][x] = -1;
        norm_vecs[3][y] = -1;

        for bund in flowpipe:
            bund_sys = bund.getIntersect()

            supp_points = np.asarray([bund_sys.max_opt(vec).x for vec in norm_vecs])
            x_points = supp_points[:, x]
            y_points = supp_points[:, y]

            ax.scatter(x_points, y_points, label=flowpipe.label, color="r")

            'Add to trajectories.'
            for p_idx, (x1, y1) in enumerate(zip(x_points, y_points)):
                supp_traj_points[p_idx][0].append(x1)
                supp_traj_points[p_idx][1].append(y1)

        'Plot the curves showing phase trajectories'
        for points in supp_traj_points:
            ax.plot(points[0], points[1], color=f"C{flow_idx}")


"""
TODO document expected dict setup.
"""


class VolumeSubplot(Subplot):

    def __init__(self, model, flowpipes, num_steps, accum_flag, plot_all_vol):
        self.accum_flag = accum_flag
        self.plot_all_vol = plot_all_vol  # Flag for plotting volumes for all methods (ConvHull, Enveloping box)
        super().__init__(model, "Volume", flowpipes, 1, num_steps)

    """
    Plots volume estimation data from self.flowpipes into input Axis object
    @params ax: Axis object to plot volume data into
    """

    def plot(self, ax):
        assert len(ax) == 1, "Only one axis object for plotting phase."
        ax = ax[0]

        axis_patches = []
        for flow_idx, flowpipe in enumerate(self.flowpipes):
            patches = self.__plot_volume_and_legend(ax, flowpipe, flow_idx)
            axis_patches += patches

        ax.set_xlabel("Time steps")
        ax.set_ylabel("Volume")
        ax.set_title(f"Volume Plot for {self.model.name}")

        ax.legend(handles=axis_patches)

    def __plot_volume_and_legend(self, ax, flowpipe, flow_idx):
        t = np.arange(0, self.num_steps, 1)
        vol_data = flowpipe.get_volume_data(accum=self.accum_flag)

        tot_vol_tup = flowpipe.all_total_volume
        tot_conv_hull_vol = tot_vol_tup.TotalFlowpipeConvHullVol
        tot_envelop_box_vol = tot_vol_tup.TotalFlowpipeEnvelopBoxVol

        if self.plot_all_vol:
            ax.plot(t, vol_data.FlowpipeConvHullVol, color=f"C{flow_idx}")
            ax.plot(t, vol_data.FlowpipeEnvelopBoxVol, color=f"C{flow_idx + self.num_flowpipes}")

            return [pat.Patch(color=f"C{flow_idx}",
                              label=f"{flowpipe.label} (Convex)"),
                    pat.Patch(color=f"C{flow_idx + self.num_flowpipes}",
                              label=f"{flowpipe.label} (EnvelopBox)")]

        else:
            if tot_conv_hull_vol > 0:
                ax.plot(t, vol_data.FlowpipeConvHullVol, color=f"C{flow_idx}")
                return [pat.Patch(color=f"C{flow_idx}",
                                  label=f"{flowpipe.label} (Convex)")]

            ax.plot(t, vol_data.FlowpipeEnvelopBoxVol, color=f"C{flow_idx}")

            return [pat.Patch(color=f"C{flow_idx}",
                              label=f"{flowpipe.label} (EnvelopBox)")]


class InitVolReachVolPlot(Subplot):

    def __init__(self, model, flowpipes, num_steps, flowpipe_indepen_data, log_scale_flag):
        super().__init__(model, "InitVolReachVol", flowpipes, 1, num_steps)
        self.flowpipe_indepen_data = flowpipe_indepen_data
        self.log_scale_flag = log_scale_flag

    def plot(self, ax):
        assert len(ax) == 1, "Only one axis object for plotting InitVolReachVol."
        ax = ax[0]

        ax.set_xlabel("Volume of Initial Box", fontsize=PlotSettings.PlotFont)
        ax.set_ylabel("Total Volume of Reachable Set", fontsize=PlotSettings.PlotFont)
        ax.set_title(f"Initial Box Volume VS Reachable Set Volume Plot for {self.model.name}",
                     fontsize=PlotSettings.PlotFont)

        ax.tick_params(axis='both', labelsize=PlotSettings.PlotFont)

        key = lambda pipe: pipe.label
        flow_dict = groupby(sorted(self.flowpipes, key=key), key=key)

        axis_patches = []
        for label_idx, (label, flowpipes_grouper) in enumerate(flow_dict):
            flowpipes = list(flowpipes_grouper)  # Grouper needs to be turned into a list
            tot_labeled_flow = len(flowpipes)

            init_vol_arr = np.empty(tot_labeled_flow)
            reach_vol_arr = np.empty(tot_labeled_flow)

            for flow_idx, flowpipe in enumerate(flowpipes):
                init_vol_arr[flow_idx] = flowpipe.init_box_volume
                reach_vol_arr[flow_idx] = flowpipe.total_volume

            ax.plot(init_vol_arr, reach_vol_arr, marker='o', color=f"C{label_idx}")
            axis_patches.append(pat.Patch(color=f"C{label_idx}",
                                          label=label))

        if self.flowpipe_indepen_data:
            label, init_vol_arr, reach_vol_arr = zip(*self.flowpipe_indepen_data)
            ax.plot(init_vol_arr, reach_vol_arr, marker='o', color=f"C{2}")
            axis_patches.append(pat.Patch(color=f"C{2}",
                                          label=label[0]))

        ax.legend(handles=axis_patches)
        y_start, y_end = ax.get_ylim()

        x_start = max(init_vol_arr[0], 0)
        x_end = init_vol_arr[-1]

        y_start = max(y_start, 0)

        num_steps = PlotSettings.NumSteps
        x_step_size = (x_end - x_start) / num_steps
        y_step_size = (y_end - y_start) / num_steps

        x_ticks = np.arange(x_start, x_end + x_step_size, x_step_size)
        y_ticks = np.arange(y_start, y_end + y_step_size, y_step_size)

        ax.set_xticks(x_ticks)
        ax.set_yticks(np.round(y_ticks))

        if self.log_scale_flag:
            ax.set_yscale('log')


"""
Object containing matplotlib figure and relevant settings and data along one axis.
"""


class Plot:

    def __init__(self):
        self.flowpipes = []
        self.trajs = []
        self.model = None
        self.num_steps = -1

    """
    Add plottable object.
    @params plottable: Traj or Flowpipe object to plot.
    """

    def add(self, plottable):
        if isinstance(plottable, TrajCollection):
            self.__add_traj(plottable)
        elif isinstance(plottable, FlowPipe):
            self.__add_flowpipe(plottable)
        else:
            raise RuntimeError("Object is not of a plottable type.")

    def plot(self, *subplots):
        axes = []
        subplot_objs = []
        for subplot_dict in subplots:
            subplot_type = subplot_dict['type']

            if subplot_type == "Projection":
                subplot = ProjectionSubplot(self.model,
                                            subplot_dict['vars'],
                                            self.flowpipes,
                                            self.num_steps,
                                            self.trajs,
                                            subplot_dict['plot_width_flag'],
                                            subplot_dict['xlims'],
                                            subplot_dict['ylims'],
                                            None if 'scale_factor' not in subplot_dict else subplot_dict['scale_factor'])

            elif subplot_type == "Phase":
                subplot = PhaseSubplot(self.model,
                                       subplot_dict['vars'],
                                       self.flowpipes,
                                       self.num_steps,
                                       subplot_dict['separate_flag'],
                                       self.trajs,
                                       subplot_dict['xlims'],
                                       subplot_dict['ylims'])

            elif subplot_type == "Volume":
                subplot = VolumeSubplot(self.model,
                                        self.flowpipes,
                                        self.num_steps,
                                        subplot_dict['accum_flag'],
                                        subplot_dict['plot_all_vol_flag'])

            elif subplot_type == "InitVolReachVol":
                subplot = InitVolReachVolPlot(self.model,
                                              self.flowpipes,
                                              self.num_steps,
                                              subplot_dict['flowpipe_indepen_data'],
                                              subplot_dict['log_scale_flag'])

            elif subplot_type == "CovidProj":
                subplot = CovidProjectionDateSubplot(self.model,
                                                     self.flowpipes,
                                                     self.num_steps,
                                                     self.trajs,
                                                     subplot_dict['data_dict'],
                                                     subplot_dict['steps_in_day'],
                                                     subplot_dict['total_pop'],
                                                     subplot_dict['time_interval'],
                                                     subplot_dict['category'])
            else:
                raise RuntimeError("Subplot type string not valid.")

            subplot_objs.append(subplot)

        figure = plt.figure(figsize=PlotSettings.fig_size, dpi=100)
        total_num_subplots = sum(subplot.num_plots for subplot in subplot_objs)

        for subplot_idx, _ in enumerate(subplot_objs):
            offset = sum(subplot.num_plots for subplot in subplots[:subplot_idx])
            subplot_tup = [figure.add_subplot(1, total_num_subplots, offset + (i + 1)) for i in
                           range(subplot.num_plots)]
            axes.append(subplot_tup)

        for axes_tup, subplot_obj in zip(axes, subplot_objs):
            subplot_obj.plot(axes_tup)

        # figure_name = "Kaa{}Phase{}--{}.png".format(self.model.name, self.__create_strat_str())

        self.__plot_figure(figure, "Fig")

    """
    Adds trajectory to be plotted.
    @params traj: Traj object to plot.
    """

    def __add_traj(self, traj_col):
        assert isinstance(traj_col,
                          TrajCollection), "Only TrajCollection objects can be added through Plot.__add_flowpipe"

        self.trajs = traj_col
        self.model = traj_col.model if self.model is None else self.model
        self.num_steps = max(self.num_steps, traj_col.traj_len)

    """
    Adds flowpipe to be plotted.
    @params flowpipe: Flowpipe object to plot.
            label: optional string argument to label the Flowpipe object in the matplotlib figure.
    """

    def __add_flowpipe(self, flowpipe):
        assert isinstance(flowpipe, FlowPipe), "Only FlowPipe objects can be added through Plot.__add_flowpipe"

        self.flowpipes.append(flowpipe)
        self.model = flowpipe.model if self.model is None else self.model
        self.num_steps = max(self.num_steps, len(flowpipe))

    """
    Creates simple string of strategy names to denote each flowpipe in an experiment. Is/was used to create unique filenames for experiment results.
    @returns string
    """

    def __create_strat_str(self):
        return ' vs '.join([str(pipe) for pipe in self.flowpipes])

    """
    Creates simple string of variable names uppercased and concatenated.
    @params var_idxs: list or tuple of variable indices pointing to its location in model.vars
    @returns string
    """

    def __create_var_str(self, var_idxs):
        return ''.join([str(self.model.vars[var_idx]).upper() for var_idx in var_idxs])

    """
    Saves or plots existing figure
    @params figure: Figure object carrying data to plot/save
            filename: filename string to save to disk
    """

    def __plot_figure(self, figure, filename):
        fig_path = self.__gen_plot_directory()
        if PlotSettings.save_fig:
            figure.savefig(os.path.join(fig_path, filename), format='png')
        else:
            plt.show()

    """
    Generates directory path used to save spreadsheet into disk.
    @returns total path to data directory
    """

    def __gen_plot_directory(self):
        data_pwd = os.path.join(PlotSettings.default_fig_path, 'Plots', date.today().isoformat(), str(self.model))
        Path(data_pwd).mkdir(parents=True, exist_ok=True)
        return data_pwd


class CombinedPlot:

    def __init__(self):
        self.intersectPlot = None
        self.sampledTrajPlot = None


class TrajLine:

    def __init__(self, line, start_point, end_point):
        self.line = line
        self.start_point = start_point
        self.end_point = end_point

    def remove(self):
        self.line.remove()
        self.start_point.remove()
        self.end_point.remove()


class SlideCompareAnimation(Plot):

    def __init__(self, *flowpipes):
        super().__init__()
        self.flowpipes = flowpipes
        self.model = flowpipes[0].model

    def animate(self, x, y, ptope_order, plot_samp_pts_flags, filename, show_recent=True):
        assert len(plot_samp_pts_flags) == len(
            self.flowpipes), "There should be a plot points Boolean flag for each flowpipe to be plotted."

        figure, ax_list = self.__init_subplots(x, y, plot_samp_pts_flags)

        'Fetch all the trajectory data for plotting.'
        num_steps = min([len(flowpipe) for flowpipe in self.flowpipes]) - 1

        'Gets the ith bundle and fetches the ptope associated to the supplied ptope order input'
        lifespan = max([flowpipe.strat.lifespan for flowpipe in self.flowpipes])

        prev_line_objs = None

        def update(i):
            ptope_idx = (i % lifespan) if show_recent else ptope_order

            nonlocal prev_line_objs

            gen_vecs_list = []
            line_objs_by_ax = []

            for ax_idx, ax in enumerate(ax_list):

                flowpipe = self.flowpipes[ax_idx]
                flowpipe_ptope = flowpipe[i].ptope(ptope_idx)

                flowpipe_traj_data = self.flowpipes[ax_idx].traj_data
                plot_samp_pt_flag = plot_samp_pts_flags[ax_idx]

                # if not flowpipe_ptope:
                #    continue

                ax_ptope, comple_ax_ptope = flowpipe_ptope

                if plot_samp_pt_flag:  # Check plotting flags

                    intersectPlot = ax.intersectPlot
                    sampledTrajPlot = ax.sampledTrajPlot

                    'Remove trajectory lines from last iteration.'
                    if i and prev_line_objs:
                        for line in prev_line_objs[ax_idx]: line.remove()

                    'If valid ptope index, plot the plot in blue and intersection of the other ptopes in green.'
                    if ptope_order < 0 and not show_recent:
                        total_intersect = flowpipe[i].getIntersect()
                        self.plot_halfspace(x, y, intersectPlot, total_intersect, idx_offset=0)
                    else:
                        'Plot ptope and its complement.'
                        self.plot_halfspace(x, y, intersectPlot, ax_ptope, idx_offset=0, alpha=0.7)
                        self.plot_halfspace(x, y, intersectPlot, comple_ax_ptope, idx_offset=2, alpha=0.4)

                    'Plot ptope and its sampled trajectory data.'
                    initial_points = flowpipe_traj_data.initial_points[i]
                    image_points = flowpipe_traj_data.image_points[i]

                    self.plot_halfspace(x, y, sampledTrajPlot, ax_ptope, idx_offset=0)  # Plot desired ptope.

                    line_obj_list = self.__plot_samp_trajs(sampledTrajPlot, x, y, initial_points, image_points)
                    line_objs_by_ax.append(line_obj_list)

                    self.plot_trajs(x, y, intersectPlot, num_steps=i + 1)

                else:
                    self.plot_halfspace(x, y, ax, ax_ptope, idx_offset=0, alpha=0.7)
                    self.plot_halfspace(x, y, ax, comple_ax_ptope, idx_offset=2, alpha=0.4)
                    line_objs_by_ax.append([])  # Placeholder for plots not having trajectory data plotted

                'Matrix of rows representing generator vectors for ax_ptope'
                gen_vecs_list.append(ax_ptope.generatorVecs)

            # self.__draw_comp_stats(ax_list[0], gen_vecs_list)
            prev_line_objs = line_objs_by_ax  # Store current line objects for removal during next step

        ani = animate.FuncAnimation(figure, update, frames=num_steps)

        Writer = animate.writers['ffmpeg']
        writer = Writer(fps=7, bitrate=-1)

        ani.save(os.path.join(PlotSettings.default_fig_path, filename + ".mp4"), writer=writer)  # save animation

    def __init_subplots(self, x, y, plot_samp_pts_flags):
        x_var, y_var = self.model.vars[x], self.model.vars[y]

        figure = plt.figure(figsize=PlotSettings.fig_size)
        num_plots = sum([1 if plot_flag else 2 for plot_flag in plot_samp_pts_flags])

        ax_list = []

        for plot_idx, plot_flag in enumerate(plot_samp_pts_flags):

            if plot_flag:
                subplot = CombinedPlot()
                subplot.intersectPlot = figure.add_subplot(1, num_plots, plot_idx + 1)
                subplot.sampledTrajPlot = figure.add_subplot(1, num_plots, plot_idx + 2)
            else:
                subplot = figure.add_subplot(1, num_plots, plot_idx + 1)

            ax_list.append(subplot)

        for ax_idx, (ax, flowpipe) in enumerate(zip(ax_list, self.flowpipes)):

            plot_flag = plot_samp_pts_flags[ax_idx]

            if plot_flag:
                ax.intersectPlot.set_xlabel(f"{x_var}")
                ax.intersectPlot.set_ylabel(f"{y_var}")
                ax.intersectPlot.set_title("Phase Plot for {}".format(flowpipe.label))

                ax.sampledTrajPlot.set_xlabel(f"{x_var}")
                ax.sampledTrajPlot.set_ylabel(f"{y_var}")
                ax.sampledTrajPlot.set_title("Phase Plot for {}".format(flowpipe.label))

            else:
                ax.set_xlabel(f"{x_var}")
                ax.set_ylabel(f"{y_var}")
                ax.set_title("Phase Plot for {}".format(flowpipe.label))

        return figure, ax_list

    def __plot_samp_trajs(self, ax, x, y, init_pts, img_pts):
        'Plot sampled trajectory lines'
        line_obj_list = []
        for init_pt, img_pt in zip(init_pts, img_pts):
            traj_line_pts = np.asarray([init_pt, img_pt])

            line, = ax.plot(traj_line_pts[:, x], traj_line_pts[:, y], c="r", linewidth=0.3)  # Plot lines.
            start_pt = ax.scatter(traj_line_pts[0, x], traj_line_pts[0, y])
            end_pt = ax.scatter(traj_line_pts[1, x], traj_line_pts[1, y], marker='x')

            traj_line_obj = TrajLine(line,
                                     start_pt,
                                     end_pt)

            line_obj_list.append(traj_line_obj)

        return line_obj_list

    def __draw_comp_stats(self, ax, gen_vecs_list):
        norm_val = []

        for vec1, vec2 in zip(*gen_vecs_list):
            vec_diff = np.subtract(vec1, vec2)
            norm_val.append(np.linalg.norm(vec_diff))

        max_diff = max(norm_val)

        patch = [pat.Patch(color='g', label=f"Largest Deviaton: {max_diff}")]
        ax.legend(handles=patch, loc='upper right', bbox_to_anchor=(1.15, 1))


class TempAnimation(Plot):

    def __init__(self, flowpipe):
        self.flowpipe = flowpipe
        self.model = flowpipe.model

    def animate(self, x, y, *strats):
        x_var, y_var = self.model.vars[x], self.model.vars[y]

        figure = plt.figure(figsize=PlotSettings.fig_size)
        ax = figure.add_subplot(1, 1, 1)

        ax.set_xlabel(f"{x_var}")
        ax.set_ylabel(f"{y_var}")
        ax.set_title("Phase Plot for {}".format(self.model.name))
        self.__draw_animation_legend(ax, *strats)

        strat_ptope_list = list(zip(*[self.flowpipe.get_strat_flowpipe(strat) for strat in strats]))

        def update(i):
            if None not in strat_ptope_list[i]:
                for ptope_idx, ptope in enumerate(strat_ptope_list[i]):
                    self.plot_halfspace(x, y, ax, ptope, idx_offset=ptope_idx)

        ani = animate.FuncAnimation(figure, update, frames=len(self.flowpipe))

        Writer = animate.writers['ffmpeg']
        writer = Writer(fps=7, bitrate=-1)

        filename = f"{self.model.name}: {' vs '.join(map(str, strats))}"
        ani.save(os.path.join(PlotSettings.default_fig_path, filename + ".mp4"), writer=writer)  # save animation

    def __draw_animation_legend(self, ax, *strats):
        axis_patches = []
        for strat_idx, strat in enumerate(strats):
            axis_patches.append(pat.Patch(color=f"C{strat_idx}", label=str(strat)))

        ax.legend(handles=axis_patches)
