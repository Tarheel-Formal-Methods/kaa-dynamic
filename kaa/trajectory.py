import numpy as np

# Create plot namespace with plotting capabilites for flowpipe/traj
# For input (flowpipes, traj) - plot2D function which takes flwopipe projections and traj and plots all of them

"""
Wrapper around list for representing arbitrary trajectories of a system.
"""
class Traj:

    def __init__(self, model):

        self.model = model
        self.vars = model.vars
        self.traj_set = {}
        self.num_points = 0

        for var in self.vars:
            self.traj_set[var] = []
    """
    Add a point from the system to the trajectory.
    @params traj_point: point to add to the trajectory.
    """
    def add_point(self, traj_point):

        assert len(traj_point) == len(self.vars), "Trajectory dimensions should match system dimensions."
        
        for var_ind, var in enumerate(self.vars):
            self.traj_set[var].append(traj_point[var_ind])
        self.num_points += 1

    """
    Returns the projection of the trajectory onto an variable axis
    @params var: var to project onto.
    @returns projection onto axis determined by var.
    """
    def get_proj(self, var):
        return self.traj_set[var]

    def get_mat(self):
        dim = len(self.vars)
        mat = np.empty((self.num_points, dim))

        for i in range(self.num_points):
            mat[i] = [ self.traj_set[var][i] for var in self.vars  ]

    def __getitem__(self, index):
        return [ (self.traj_set[var])[index] for var in self.vars ]

    def __len__(self):
        return self.num_points

    @property
    def model_name(self):
        return self.model.name
