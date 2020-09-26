import matplotlib as plt
import numpy as np
from scipy.spatial import HalfspaceIntersection
import os

from kaa.settings import PlotSettings
from kaa.trajectory import Traj
from kaa.flowpipe import FlowPipe
from kaa.timer import Timer
from kaa.lputil import minLinProg, maxLinProg


plt.rcParams.update({'font.size': PlotSettings.plot_font})

"""
Object containing matplotlib figure and relevant settings and data along one axis.
"""
class Plot:

    def __init__(self):

        self.flowpipes = []
        self.trajs = []
        self.model = None
        self.num_steps = 0
    """
    Add plottable object.
    @params plottable: Traj or Flowpipe object to plot.
    """
    def add(self, plottable, label=None):

        if isinstance(plottable, Traj):
            self.__add_traj(plottable)
        elif isinstance(plottable, FlowPipe):
            self.__add_flowpipe(plottable, label=label)
        else:
            raise RuntimeError("Object is not of a plottable type.")
        
    """
    Adds trajectory to be plotted.
    @params traj: Traj object to plot.
    """
    def __add_traj(self, traj):

        assert isinstance(traj, Traj), "Only Traj objects can be added through Plot.__add_flowpipe"
        if self.model is not None:
            assert self.model.name == traj.model_name, "Trajectories and Plot must describe the same system."

        self.trajs.append(traj)
        self.model = traj.model if self.model is None else self.model
        self.num_steps = max(self.num_steps, len(traj))
        
    """
    Adds flowpipe to be plotted.
    @params flowpipe: Flowpipe object to plot.
            label: optional string argument to label the Flowpipe object in the matplotlib figure.
    """
    def __add_flowpipe(self, flowpipe, label=None):

        assert isinstance(flowpipe, FlowPipe), "Only FlowPipe objects can be added through Plot.__add_flowpipe"
        if self.model is not None:
              assert self.model.name == flowpipe.model_name, "FlowPipe and Plot must describe the same system."

        self.flowpipes.append((label, flowpipe))
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

        figure = plt.figure.Figure(figsize=PlotSettings.fig_size)
        figure.subplots_adjust(hspace=0.3,wspace=0.2)

        'Hackish way of adding subplots to Figure objects.'
        ax = [ figure.add_subplot(1, num_var, i+1) for i in range(num_var) ] if overlap else \
             [[figure.add_subplot(num_flowpipes, num_var, (y*num_var)+(i+1)) for i in range(num_var)] for y in range(num_flowpipes)]

        for ax_idx, var_ind in enumerate(var_tup):

            var = self.model.vars[var_ind]
            name = self.model.name
            t = np.arange(0, self.num_steps, 1)

            for traj_idx, traj in enumerate(self.trajs):
                x = np.arange(len(traj))
                y = traj.get_traj_proj(var)
                ax[ax_idx].plot(x, y, color="C{}".format(traj_idx))

            for flow_idx, (label, flowpipe) in enumerate(self.flowpipes):
                flow_min, flow_max = flowpipe.get2DProj(var_ind)
                flowpipe_label = name if label is None else label

                curr_ax = ax[ax_idx] if overlap else ax[flow_idx][ax_idx]

                curr_ax.fill_between(t, flow_min, flow_max, label=flowpipe_label, color="C{}".format(flow_idx), alpha=0.5)
                curr_ax.set_xlabel("t: time steps")
                curr_ax.set_ylabel(("Reachable Set for {}".format(var)))
                curr_ax.set_title("Projection of Reachable Set for {} Variable: {}".format(name, var))
                curr_ax.legend()

        if PlotSettings.save_fig:
            var_str = ''.join([str(self.model.vars[var_idx]).upper() for var_idx in var_tup])
            strat_str = ' vs '.join([str(pipe) for label, pipe in self.flowpipes])
            
            figure_name = "Kaa{}Proj{}--{}.png".format(self.model.name, var_str, strat_str)
            figure.savefig(os.path.join(path, figure_name), format='png')
        else:
            figure.show()

    """

    Plots phase between two variables of dynamical system.
    NOTE: This method for now creates rather crude, segmented phase plots by simply scattering the support points.

    @params x: index of variable to be plotted as x-axis of desired phase
            y: index of variable to be plotted as y-axis of desired phase
    """
    def plot2DPhase(self, x, y):

        Timer.start('Phase')

        dim = self.model.dim
        x_var, y_var = self.model.vars[x], self.model.vars[y]

        'Define the following projected normal vectors.'
        norm_vecs = np.zeros([4, dim])
        norm_vecs[0][x] = 1; norm_vecs[1][y] = 1;
        norm_vecs[2][x] = -1; norm_vecs[3][y] = -1;
    
        figure = plt.figure.Figure(figsize=PlotSettings.fig_size)
        ax = figure.add_subplot(1,1,1)

        comple_dim = [i for i in range(dim) if i not in [x,y]]

        'Initialize objective function for Chebyshev intersection LP routine.'
        c = [0 for _ in range(dim + 1)]
        c[-1] = 1

        for flow_idx, (flow_label, flowpipe) in enumerate(self.flowpipes):

            for bund in flowpipe:
                bund_A, bund_b = bund.getIntersect()

                supp_points = np.asarray([ maxLinProg(vec, bund_A, bund_b).x  for vec in norm_vecs ])
                x_points = supp_points[:,x]
                y_points = supp_points[:,y]

                ax.scatter(x_points, y_points, label=flow_label, color="C{}".format(flow_idx))

            ax.set_xlabel(f'{x_var}')
            ax.set_ylabel(f'{y_var}')
            ax.set_title("Projection of Phase Plot for {} Variables: {}".format(self.model.name, (x_var, y_var)))

        if PlotSettings.save_fig:
            var_str = ''.join([str(self.model.vars[var]).upper() for var in [x,y]])
            figure_name = "Kaa{}Phase{}.png".format(flowpipe.model.name, var_str)

            figure.savefig(os.path.join(PlotSettings.default_fig_path, figure_name), format='png')
        else:
            figure.show()

        phase_time = Timer.stop('Phase')
        print("Plotting phase for dimensions {}, {} done -- Time Spent: {}".format(x_var, y_var, phase_time))


    """
    Use scipy.HalfspaceIntersection to fill in phase plot projections.
    """
    def __halfspace_inter_plot(bund, x, y, ax, norm_vecs):
        'Compute the normal vector offsets'
        bund_off = np.empty([len(norm_vecs), 1])

        for i in range(len(norm_vecs)):
            bund_off[i] = minLinProg(np.negative(norm_vecs[i]), bund_A, bund_b).fun

        'Remove irrelevant model.dimensions. Mostly doing this to make HalfspaceIntersection happy.'
        phase_intersect = np.hstack((norm_vecs, bund_off))
        phase_intersect = np.delete(phase_intersect, comple_dim, axis=1)

        'Compute Chebyshev center of intersection.'
        row_norm = np.reshape(np.linalg.norm(norm_vecs, axis=1), (norm_vecs.shape[0], 1))
        center_A = np.hstack((norm_vecs, row_norm))

        neg_bund_off = np.negative(bund_off)
        center_pt = maxLinProg(c, center_A, list(neg_bund_off.flat)).x
        center_pt = np.asarray([b for b_i, b in enumerate(center_pt) if b_i in [x, y]])

        'Run scipy.spatial.HalfspaceIntersection.'
        hs = HalfspaceIntersection(phase_intersect, center_pt)
        inter_x, inter_y = zip(*hs.intersections)
        ax.fill(inter_x, inter_y, label=flow_label, color="C{}".format(flow_idx))
