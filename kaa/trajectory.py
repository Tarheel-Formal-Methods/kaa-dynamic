import numpy as np

"""
Wrapper around list for representing arbitrary trajectories of a system.
"""
class Traj:

    def __init__(self, model, initial_point, steps=0):

        self.model = model
        self.vars = model.vars
        self.traj_set = {}
        self.num_points = 0

        for var in self.vars:
            self.traj_set[var] = []

        'Initialize point and propagate for steps.'
        self.add_point(initial_point)
        self.propagate(steps)
        
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
    Propagate tip of trajectory for alloted number of time_steps
    @params time_steps: number of time steps to generate trajectory
    """
    def propagate(self, time_steps):

        df = self.model.f

        'Propagate the points according to the dynamics for designated number of time steps.'
        prev_point = self.end_point

        for _ in range(time_steps):
            
            var_sub = [ (var, prev_point[var_idx]) for var_idx, var in enumerate(self.vars)]
            next_point = [ f.subs(var_sub) for f in df ]
            self.add_point(next_point)

            prev_point = next_point

    """
    Returns the projection of the trajectory onto an variable axis
    @params var: var to project onto.
    @returns projection onto axis determined by var.
    """
    def get_proj(self, var):
        return self.traj_set[var]

    """
    Returns numpy matrix with rows containing trajectory points.
    @returns matrix containing trajectory points.
    """
    def get_mat(self):
        dim = len(self.vars)
        mat = np.empty((self.num_points, dim))

        for i in range(self.num_points):
            mat[i] = [ self.traj_set[var][i] for var in self.vars  ]

    @property
    def end_point(self):

        assert self.num_points != 0, "Trajectory must be popoulated before querying its terminal point."
        return self[-1]

    @property
    def start_point(self):
        assert self.num_points != 0, "Trajectory must be popoulated before querying its terminal point."
        return self[0]

    @property
    def model_name(self):
        return self.model.name
    

    def __getitem__(self, index):
        return [ (self.traj_set[var])[index] for var in self.vars ]

    def __len__(self):
        return self.num_points
