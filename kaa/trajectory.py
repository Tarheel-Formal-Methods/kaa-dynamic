import matplotlib as plt


# Create plot namespace with plotting capabilites for flowpipe/traj
# For input (flowpipes, traj) - plot2D function which takes flwopipe projections and traj and plots all of them

class 2DTrajectory:

    def __init__(self, vars):

        self.vars = vars
        self.traj_set = {}

        for var in vars:
            self.traj_set[var] = set()

    def add_traj(self, var, traj):
        traj_set[var].add(traj)


    def plotTraj(self):
