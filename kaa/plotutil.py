import matplotlib as plt
import numpy as np
import os

from kaa.settings import PlotSettings
from kaa.trajectory import Traj
from kaa.flowpipe import FlowPipe

"""
Object containing matplotlib figure and relevant settings and data along one axis.
"""

class Plot:

    def __init__(self):

        self.figure = plt.figure.Figure(figsize = PlotSettings.fig_size)
        self.ax = self.figure.add_subplot(111)
        
        self.flowpipes = []
        self.trajs = []
        self.model = None
        self.num_steps = 0

    def add(self, plottable):

        if isinstance(plottable, Traj):
            self.__add_traj(plottable)
        elif isinstance(plottable, FlowPipe):
            self.__add_flowpipe(plottable)

    def __add_traj(self, traj_data):

        #assert isinstance(traj_data, Traj), "Only Traj objects can be added through Plot.add_flowpipe"
        if self.model is not None:
            assert isinstance(self.model, traj_data.model), "Trajectories and Plot must describe the same system."

        self.trajs.append(traj_data)
        self.model = traj_data.model if self.model is None else self.model
        self.num_steps = max(self.num_steps, len(traj_data))

    def __add_flowpipe(self, flowpipe):

        #assert isinstance(flowpipe, FlowPipe), "Only FlowPipe objects can be added through Plot.add_flowpipe"
        if self.model is not None:
              assert isinstance(self.model, flowpipe.model), "FlowPipe and Plot must describe the same system."
        
        self.flowpipes.append(flowpipe)
        self.model = flowpipe.model if self.model is None else self.model
        self.num_steps = max(self.num_steps, len(flowpipe))

    def plot(self, var_ind):

        assert self.model is not None, "No data has been added to Plot."

        var = self.model.vars[var_ind]
        name = self.model.name
        t = np.arange(0, self.num_steps, 1)

        for traj in self.trajs:
            x = np.arange(len(traj))
            y = traj.get_traj_proj(var)
            self.ax.plot(x, y, color="C{}".format(len(self.trajs)-1))

        for flowpipe in self.flowpipes:
            flow_min, flow_max = flowpipe.get2DProj(var_ind)
            self.ax.fill_between(t, flow_min, flow_max, color="C{}".format(len(self.flowpipes)-1))

        self.ax.set_xlabel("t: time steps")
        self.ax.set_ylabel(("Reachable Set for {}".format(var)))
        self.ax.set_title("Projection of Reachable Set for {} Variable: {}".format(name, var))

        if PlotSettings.save_fig:
            var_str = ''.join([str(var).upper()])
            figure_name = "Kaa{}Proj{}.png".format(self.model.name, var_str)

            self.figure.savefig(os.path.join(PlotSettings.fig_path, figure_name), format='png')
        else:
            self.figure.show()
