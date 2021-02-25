import os
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as pat
import matplotlib.animation as animate
from matplotlib.legend import Legend

from kaa.settings import PlotSettings
from kaa.trajectory import TrajCollection, Traj
from kaa.flowpipe import FlowPipe
from kaa.timer import Timer
from kaa.parallelotope import LinearSystem

plt.rcParams.update({'font.size': PlotSettings.plot_font})

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
    def add(self, plottable, label=None):
        if isinstance(plottable, TrajCollection):
            self.__add_traj(plottable)
        elif isinstance(plottable, FlowPipe):
            self.__add_flowpipe(plottable)
        else:
            raise RuntimeError("Object is not of a plottable type.")

    """
    Adds trajectory to be plotted.
    @params traj: Traj object to plot.
    """
    def __add_traj(self, traj_col):
        assert isinstance(traj_col, TrajCollection), "Only TrajCollection objects can be added through Plot.__add_flowpipe"
        #if self.model is not None:
        #    assert self.model.name == traj.model_name, "Trajectories and Plot must describe the same system."

        self.trajs = traj_col
        self.model = traj_col.model if self.model is None else self.model
        self.num_steps = max(self.num_steps, traj_col.max_traj_len)

    """
    Adds flowpipe to be plotted.
    @params flowpipe: Flowpipe object to plot.
            label: optional string argument to label the Flowpipe object in the matplotlib figure.
    """
    def __add_flowpipe(self, flowpipe):

        assert isinstance(flowpipe, FlowPipe), "Only FlowPipe objects can be added through Plot.__add_flowpipe"
        #if self.model is not None:
        #      assert self.model.name == flowpipe.model_name, "FlowPipe and Plot must describe the same system."

        self.flowpipes.append(flowpipe)
        self.model = flowpipe.model if self.model is None else self.model
        self.num_steps = max(self.num_steps, len(flowpipe))

    """
    Plots the projections of the trajectories and flowpipes stored in Plot object.
    @params: *var_tup: indices of desired variables
             path: optional variable designated the path to store the generated matplotlib figure.
    """
    def plot(self, *var_tup, path=PlotSettings.default_fig_path, overlap=True):
        assert self.model is not None, "No data has been added to the Plot object."
        num_var = len(var_tup)
        num_flowpipes = len(self.flowpipes)

        figure = plt.figure(figsize=PlotSettings.fig_size)
        figure.subplots_adjust(hspace=0.3,wspace=0.2)

        'Hackish way of adding subplots to Figure objects.'
        ax = [ figure.add_subplot(1, num_var, i+1) for i in range(num_var) ] if overlap else \
             [[figure.add_subplot(num_flowpipes, num_var, (y*num_var)+(i+1)) for i in range(num_var)] for y in range(num_flowpipes)]

        for ax_idx, var_ind in enumerate(var_tup):

            var = self.model.vars[var_ind]
            name = self.model.name
            t = np.arange(0, self.num_steps, 1)

            self.__plot_trajs(ax[ax_idx]) #fix this

            for flow_idx, flowpipe in enumerate(self.flowpipes):
                flow_min, flow_max = flowpipe.get2DProj(var_ind)
                flowpipe_label = flowpipe.label

                curr_ax = ax[ax_idx] if overlap else ax[flow_idx][ax_idx]
                curr_ax.fill_between(t, flow_min, flow_max, label=flowpipe_label, color="C{}".format(flow_idx), alpha=0.5)

                curr_ax.set_xlabel("t: time steps")
                curr_ax.set_ylabel(f"Reachable Set for {var}")
                curr_ax.set_title(f"Projection of Reachable Set for {name} Variable: {var}")
                curr_ax.legend()

        figure_name = "Kaa{}Proj{}--{}.png".format(self.model.name, self.__create_var_str(var_tup), self.__create_strat_str())

        self.__plot_figure(figure, figure_name)

    """
    Plots phase between two variables of dynamical system.
    NOTE: This method for now creates rather crude, segmented phase plots by simply scattering the support points.

    @params x: index of variable to be plotted as x-axis of desired phase
            y: index of variable to be plotted as y-axis of desired phase
    """
    def plot2DPhase(self, x, y, separate, lims=None, volume_chart=False):
        assert len(self.flowpipes) != 0, "Plot Object must have at least one flowpipe to plot for 2DPhase."

        Timer.start('Phase')

        dim = self.model.dim

        figure = plt.figure(figsize=PlotSettings.fig_size)
        phase_ax = figure.add_subplot(1,2,1)
        vol_ax = figure.add_subplot(1,2,2)

        for flow_idx, flowpipe in enumerate(self.flowpipes):
            self.__halfspace_inter_plot(flowpipe, flow_idx, x, y, phase_ax, separate)
            #self.__support_plot(flowpipe, flow_idx,  x, y, ax)

        self.plot_trajs(x, y, phase_ax)
        self.__phase_plot_legend(x, y, phase_ax, lims)

        if volume_chart:
            'Add volume data'
            self.plot_volume(vol_ax)

        figure_name = "Kaa{}Phase{}.png".format(flowpipe.model.name, self.__create_var_str([x,y]))
        self.__plot_figure(figure, figure_name)

        phase_time = Timer.stop('Phase')
        print("Plotting phase for dimensions {}, {} done -- Time Spent: {}".format(x_var, y_var, phase_time))

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
    Routine to populate the legend information for a phase plot.
    @params x: index of x variable
            y: index of y variable
            phase_ax: Axis object associated to phase plot.
            lims: optional data denoting x,y axis dimensions
    """
    def __phase_plot_legend(self, x, y, phase_ax, lims):
        x_var, y_var = self.model.vars[x], self.model.vars[y]

        phase_ax.set_xlabel(f'{x_var}')
        phase_ax.set_ylabel(f'{y_var}')
        phase_ax.set_title("Projection of Phase Plot for {} Variables: {}".format(self.model.name, (x_var, y_var)))

        axis_patches = []
        for flow_idx, flowpipe in enumerate(self.flowpipes):
            for strat_idx, strat in enumerate(flowpipe.strats):
                axis_patches.append(pat.Patch(color = f"C{flow_idx + strat_idx}", label=str(strat)))

        phase_ax.legend(handles=axis_patches)

        if lims is not None:
           phase_ax.set_xlim(lims)
           phase_ax.set_ylim(lims)

    """
    Routine to plot trajectory data (self.trajs) into an input Axis object
    @params x: index of x variable
            y: index of y variable
            ax: Axis object to plot trajectory data into
    """
    def plot_trajs(self, x, y, ax, num_steps=0):
        if not self.trajs: return

        x_var, y_var = self.model.vars[x], self.model.vars[y]
        
        for traj_idx, traj in enumerate(self.trajs):
            x_coord = traj.get_proj(x_var)
            y_coord = traj.get_proj(y_var)

            if num_steps:
                x_coord = x_coord[:num_steps]
                y_coord = y_coord[:num_steps]

            ax.plot(x_coord, y_coord, color="k", linewidth=1)
            ax.scatter(x_coord, y_coord, color=f"C{traj_idx}", s=1)


    def support_plot(self, flowpipe, flow_idx, x, y, ax):
        dim = self.model.dim

        'Define the following projected normal vectors.'
        norm_vecs = np.zeros([4, dim])
        norm_vecs[0][x] = 1; norm_vecs[1][y] = 1;
        norm_vecs[2][x] = -1; norm_vecs[3][y] = -1;

        for bund in flowpipe:
            bund_sys = bund.getIntersect()

            supp_points = np.asarray([bund_sys.max_opt(vec).x  for vec in norm_vecs])
            x_points = supp_points[:,x]
            y_points = supp_points[:,y]

            ax.scatter(x_points, y_points, label=flowpipe.label, color="r")

            'Add to trajectories.'
            for p_idx, (x1, y1) in enumerate(zip(x_points, y_points)):
                supp_traj_points[p_idx][0].append(x1)
                supp_traj_points[p_idx][1].append(y1)

        'Plot the curves showing phase trajectories'
        for points in supp_traj_points:
            ax.plot(points[0], points[1], color=f"C{flow_idx}")


    """
    Use scipy.HalfspaceIntersection to fill in phase plot projections.
    @params: flowpipe: FlowPipe object to plot.
    """
    def __halfspace_inter_plot(self, flowpipe, flow_idx, x, y, ax, separate):
        for bund in flowpipe:
            if not separate:
                'Temp patch. Revise to start using LinearSystems for future work.'
                self.plot_halfspace(x, y, ax, bund.getIntersect(), idx_offset=flow_idx)
            else:
                for ptope_idx, ptope in enumerate(bund.all_ptopes):
                    self.plot_halfspace(x, y, ax, ptope, idx_offset=flow_idx+ptope_idx)

    """
    Plot linear system through scipy.HalfspaceIntersection
    """
    def plot_halfspace(self, x, y, ax, sys, idx_offset=0, alpha=0.4):

        'Routines in Bundle give None values for queries of non-existent ptopes'
        if sys is None: return sys

        dim = self.model.dim
        comple_dim = np.asarray([True if i in [x,y] else False for i in range(dim)])

        vertices = sys.vertices
        proj_vertices = np.unique(vertices[:,comple_dim], axis=0).tolist()
        center_pt = sys.chebyshev_center.center

        'Sort by polar coordinates to ensure proper plotting of boundary'
        proj_vertices.sort(key=lambda v: math.atan2(v[1] - center_pt[1] , v[0] - center_pt[0]))

        ptope = pat.Polygon(proj_vertices, fill=True, color=f"C{idx_offset}", alpha=alpha)
        ax.add_patch(ptope)

        inter_x, inter_y = zip(*proj_vertices)
        ax.scatter(inter_x, inter_y, s=0.1)

    """
    Plots volume estimation data from self.flowpipes into input Axis object
    @params ax: Axis object to plot volume data into
    """
    def plot_volume(self, ax):
        num_flowpipes = len(self.flowpipes)
        t = np.arange(0, self.num_steps, 1)

        axis_patches = []
        for flow_idx, flowpipe in enumerate(self.flowpipes):
            vol_data = flowpipe.get_volume_data()
            ax.plot(t, vol_data, color=f"C{flow_idx}")
            axis_patches.append(pat.Patch(color = f"C{flow_idx}", label=f"{flowpipe.label} (Total Volume: {flowpipe.total_volume})"))

        ax.set_xlabel("Time steps")
        ax.set_ylabel("Volume")
        ax.set_title("Volume Plot for {}".format(self.model.name))

        ax.legend(handles=axis_patches)

    """
    Saves or plots existing figure
    @params figure: Figure object carrying data to plot/save
            filename: filename string to save to disk
    """
    def __plot_figure(self, figure, filename):
        if PlotSettings.save_fig:
            figure.savefig(os.path.join(PlotSettings.default_fig_path, filename), format='png')
        else:
            plt.show()

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

    def animate(self, x, y, ptope_order, plot_samp_pts_flags, filename, show_recent=False):
        assert len(plot_samp_pts_flags) == len(self.flowpipes), "There should be a plot points Boolean flag for each flowpipe to be plotted."

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
                flowpipe_ptope = flowpipe[i].ptope(ptope_order)

                flowpipe_traj_data = self.flowpipes[ax_idx].traj_data
                plot_samp_pt_flag = plot_samp_pts_flags[ax_idx]

                if not flowpipe_ptope:
                    continue

                ax_ptope, comple_ax_ptope = flowpipe_ptope

                if plot_samp_pt_flag: #Check plotting flags

                    intersectPlot = ax.intersectPlot
                    sampledTrajPlot = ax.sampledTrajPlot

                    'Remove trajectory lines from last iteration.'
                    if i and prev_line_objs:
                        for line in prev_line_objs[ax_idx]: line.remove()

                    'Plot ptope and its complement.'
                    self.plot_halfspace(x, y, intersectPlot, ax_ptope, idx_offset=0)
                    self.plot_halfspace(x, y, intersectPlot, comple_ax_ptope, idx_offset=2)

                    'Plot ptope and its sampled trajectory data.'
                    initial_points = flowpipe_traj_data.initial_points[i]
                    image_points = flowpipe_traj_data.image_points[i]

                    self.plot_halfspace(x, y, sampledTrajPlot, ax_ptope, idx_offset=0) #Plot desired ptope.

                    line_obj_list = self.__plot_samp_trajs(sampledTrajPlot, x, y, initial_points, image_points)
                    line_objs_by_ax.append(line_obj_list)

                    self.plot_trajs(x, y, intersectPlot, num_steps=i+1)

                else:
                    self.plot_halfspace(x, y, ax, ax_ptope, idx_offset=0)
                    self.plot_halfspace(x, y, ax, comple_ax_ptope, idx_offset=2)
                    line_objs_by_ax.append([]) #Placeholder for plots not having trajectory data plotted

                'Matrix of rows representing generator vectors for ax_ptope'
                gen_vecs_list.append(ax_ptope.generatorVecs)

            #self.__draw_comp_stats(ax_list[0], gen_vecs_list)
            prev_line_objs = line_objs_by_ax #Store current line objects for removal during next step
                
        ani = animate.FuncAnimation(figure, update, frames=num_steps)

        Writer = animate.writers['ffmpeg']
        writer = Writer(fps=7,bitrate=-1)

        ani.save(os.path.join(PlotSettings.default_fig_path, filename + ".mp4"), writer=writer) #save animation

    def __init_subplots(self, x, y, plot_samp_pts_flags):
        x_var, y_var = self.model.vars[x], self.model.vars[y]

        figure = plt.figure(figsize=PlotSettings.fig_size)
        num_plots = sum([1 if plot_flag else 2 for plot_flag in plot_samp_pts_flags])

        ax_list = []
        
        for plot_idx, plot_flag in enumerate(plot_samp_pts_flags):

            if plot_flag:
                subplot = CombinedPlot()
                subplot.intersectPlot = figure.add_subplot(1, num_plots, plot_idx+1)
                subplot.sampledTrajPlot = figure.add_subplot(1, num_plots, plot_idx+2)
            else:
                subplot = figure.add_subplot(1, num_plots, plot_idx+1)
                
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

            line, = ax.plot(traj_line_pts[:,x], traj_line_pts[:,y], c="r") #Plot lines.
            start_pt = ax.scatter(traj_line_pts[0,x], traj_line_pts[0,y])
            end_pt = ax.scatter(traj_line_pts[1,x], traj_line_pts[1,y], marker='x', s=1.3)

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
        ax = figure.add_subplot(1,1,1)

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
        writer = Writer(fps=7,bitrate=-1)

        filename = f"{self.model.name}: {' vs '.join(map(str, strats))}"
        ani.save(os.path.join(PlotSettings.default_fig_path, filename + ".mp4"), writer=writer) #save animation

    def __draw_animation_legend(self, ax, *strats):
        axis_patches = []
        for strat_idx, strat in enumerate(strats):
            axis_patches.append(pat.Patch(color=f"C{strat_idx}", label=str(strat)))

        ax.legend(handles=axis_patches)
