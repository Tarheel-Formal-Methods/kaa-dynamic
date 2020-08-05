import matplotlib as plt
import numpy as np

from kaa.settings import PlotSettings
from kaa.trajectory import Traj
from kaa.flowpipe import FlowPipe

"""
Object containing matplotlib figure and relevant settings and data along one axis.
"""

class Plot:

    def __init__(self, model, var):

        self.figure = plt.Figure(figsize = PlotSettings.fig_size)
        self.ax = self.figure.add_subplot(111)
        self.ax.set_xlabel("t: time steps")
        self.ax.set_ylabel(("Reachable Set for {}".format(var)))
        self.ax.set_title("Projection of Reachable Set for {} Variable: {}".format(model.name, var))
        
        self.num_flowpipes = 0
        self.num_traj = 0

    def add_traj(self, traj_data):

        assert isinstance(traj_data, Traj)

        traj_list = traj_data[var]

        for traj in traj_list:
            x = np.arange(len(traj))
            self.ax.plot(x, traj, color="C{}".format(self.num_traj))
            self.num_traj += 1


    def add_flowpipe(self, flowpipe):

        assert isinstance(flowpipe, FlowPipe)

        flow_min, flow_max = flowpipe.get2DProj()
        ax.fill_between(flow_min, flow_max, color="C{}".format(self.num_flowpipes))
        self.num_flowpipes += 1

    def plot(self):

        if PlotSettings.save_fig:
            var_str = ''.join([str(self.vars[var]).upper() for var in vars_tup])
            figure_name = "Kaa{}Proj{}.png".format(self.name, var_str)

            self.figure.savefig(os.path.join(PlotSettings.fig_path, figure_name), format='png')
        else:
            self.figure.show()
