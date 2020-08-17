import matplotlib as plt

# Create plot namespace with plotting capabilites for flowpipe/traj
# For input (flowpipes, traj) - plot2D function which takes flwopipe projections and traj and plots all of them

"""
Wrapper around list for representing arbitrary trajectories of a system.
"""
class Traj:

    def __init__(self, model):

        self.vars = model.vars
        self.traj_set = {}
        self.num_traj = 0

        for var in vars:
            self.traj_set[var] = []

    def add_traj_point(self, traj_point):

        assert len(traj_point) == len(self.vars), "Trajectory dimensions should match system dimensions."
        
        for var_ind, var in self.vars:
            traj_set[var].add(traj_point[var_ind])
        self.num_traj += 1

    def get_traj_proj(self, var):
        return traj_set[var]

    def __getitem__(self, index):
        return [ (traj_set[var])[index] for var in self.vars ]

    def __len__(self):
        return self.num_traj
