import matplotlib as plt
import numpy as np
import os

from kaa.settings import PlotSettings
from kaa.trajectory import Traj
from kaa.flowpipe import FlowPipe

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

        assert isinstance(traj, Traj), "Only Traj objects can be added through Plot.add_flowpipe"
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

        assert isinstance(flowpipe, FlowPipe), "Only FlowPipe objects can be added through Plot.add_flowpipe"
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
    def plot(self, *var_tup, path=PlotSettings.default_fig_path):

        assert self.model is not None, "No data has been added to the Plot object."
        num_var = len(var_tup)
        figure = plt.figure.Figure(figsize=PlotSettings.fig_size)

        'Hackish way of adding subplots to Figure objects.'
        ax = [ figure.add_subplot(1, num_var, i+1) for i in range(num_var) ]

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

                ax[ax_idx].fill_between(t, flow_min, flow_max, label=flowpipe_label, color="C{}".format(flow_idx))

                ax[ax_idx].set_xlabel("t: time steps")
                ax[ax_idx].set_ylabel(("Reachable Set for {}".format(var)))
                ax[ax_idx].set_title("Projection of Reachable Set for {} Variable: {}".format(name, var))
                ax[ax_idx].legend()

        if PlotSettings.save_fig:
            var_str = ''.join([str(self.model.vars[var_idx]).upper() for var_idx in var_tup])
            figure_name = "Kaa{}Proj{}.png".format(self.model.name, var_str)

            figure.savefig(os.path.join(path, figure_name), format='png')
        else:
            figure.show()
