import os
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as pat
import matplotlib.animation as animate
from scipy.spatial import HalfspaceIntersection

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
    def plot2DPhase(self, x, y, separate=False, lims=None):
        assert len(self.flowpipes) != 0, "Plot Object must have at least one flowpipe to plot for 2DPhase."

        Timer.start('Phase')

        dim = self.model.dim

        figure = plt.figure(figsize=PlotSettings.fig_size)
        phase_ax = figure.add_subplot(1,2,1)
        vol_ax = figure.add_subplot(1,2,2)

        for flow_idx, flowpipe in enumerate(self.flowpipes):
            self.__halfspace_inter_plot(flowpipe, flow_idx, x, y, phase_ax, separate)
            #self.__support_plot(flowpipe, flow_idx,  x, y, ax)

        self.__plot_trajs(x, y, phase_ax)
        self.__phase_plot_legend(x, y, phase_ax, lims)

        'Add volume data'
        self.__plot_volume(vol_ax)

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
    def __plot_trajs(self, x, y, ax):
        x_var, y_var = self.model.vars[x], self.model.vars[y]

        for traj_idx, traj in enumerate(self.trajs):
            x_coord = traj.get_proj(x_var)
            y_coord = traj.get_proj(y_var)

            ax.plot(x_coord, y_coord, color=f"C{traj_idx}", linewidth=1)
            ax.scatter(x_coord, y_coord, color=f"C{traj_idx}", s=0.5)


    def __support_plot(self, flowpipe, flow_idx, x, y, ax):
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
                for ptope_idx, ptope in enumerate(bund.ptopes):
                    self.plot_halfspace(x, y, ax, ptope, idx_offset=flow_idx+ptope_idx)
                    
    """
    Plot linear system through scipy.HalfspaceIntersection
    """
    def plot_halfspace(self, x, y, ax, sys, idx_offset=0):

        dim = self.model.dim
        comple_dim = np.asarray([ True if i in [x,y] else False for i in range(dim) ])

        'Halfspace constraint matrix'
        A = sys.A
        b = sys.b

        phase_intersect = np.hstack((A, - np.asarray([b]).T))
        center_pt = np.asarray(sys.chebyshev_center.center)

        'Run scipy.spatial.HalfspaceIntersection.'
        hs = HalfspaceIntersection(phase_intersect, center_pt)
        vertices = np.asarray(hs.intersections)

        proj_vertices = np.unique(vertices[:,comple_dim], axis=0).tolist()

        'Sort by polar coordinates to ensure proper plotting of boundary'
        proj_vertices.sort(key=lambda v: math.atan2(v[1] - center_pt[1] , v[0] - center_pt[0]))

        ptope = pat.Polygon(proj_vertices, fill=True, color=f"C{idx_offset}", alpha=0.4)
        ax.add_patch(ptope)
  
        inter_x, inter_y = zip(*proj_vertices)
        ax.scatter(inter_x, inter_y, s=0.1)

    """
    Plots volume estimation data from self.flowpipes into input Axis object
    @params ax: Axis object to plot volume data into
    """
    def __plot_volume(self, ax):
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
