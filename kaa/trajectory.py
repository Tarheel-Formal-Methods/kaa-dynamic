import numpy as np
from scipy.spatial import ConvexHull

from kaa.settings import KaaSettings

from kaa.log import Output
"""
Wrapper around list for representing arbitrary trajectories of a system.
"""
class Traj:

    def __init__(self, model, initial_point, steps=0):
        self.model = model
        self.dim = model.dim
        self.vars = model.vars
        self.traj_mat = np.empty((1,self.dim))

        'Initialize point and propagate for steps.'
        self.traj_mat[0] = initial_point
        self.propagate(steps)

    @property
    def num_points(self):
        return len(self.traj_mat)

    """
    Add a point from the system to the trajectory.
    @params traj_point: point to add to the trajectory.
    """
    def add_point(self, traj_point):
        assert len(traj_point) == len(self.vars), "Trajectory dimensions should match system dimensions."
        self.traj_mat = np.append(self.traj_mat, [traj_point], axis=0)


    """
    Propagate tip of trajectory for alloted number of time_steps
    @params time_steps: number of time steps to generate trajectory
    """
    def propagate(self, time_steps):

        if KaaSettings.Parallelize:
            df = self.model.f
        else:
            df = self.model.lambdified_f

        'Propagate the points according to the dynamics for designated number of time steps.'
        prev_point = self.end_point

        for _ in range(time_steps):
            var_sub = [(var, prev_point[var_idx]) for var_idx, var in enumerate(self.vars)]

            if KaaSettings.Parallelize:
                next_point = [f.subs(var_sub) for f in df]
            else:
                next_point = df(*prev_point)
            self.add_point(next_point)

            prev_point = next_point


    """
    Returns the projection of the trajectory onto an variable axis
    @params var: var to project onto.
    @returns projection onto axis determined by var.
    """
    def get_proj(self, var_idx):
        return self.traj_mat[:,var_idx]

    @property
    def end_point(self):
        assert self.num_points != 0, "Trajectory must be popoulated before querying its terminal point."
        return self[-1]

    @property
    def start_point(self):
        assert self.num_points != 0, "Trajectory must be popoulated before querying its terminal point."
        return self[0]

    def __getitem__(self, index):
        return self.traj_mat[index]

    def __len__(self):
        return self.num_points

class TrajCollection:

    def __init__(self, model, traj_list):
        assert isinstance(traj_list, list), "input must be list of Traj objects"
        self.traj_list = traj_list
        self.num_trajs = len(self.traj_list)
        self.model = model
        self.dim = model.dim
        self.traj_len = min(traj.num_points for traj in self.traj_list)

    @property
    def start_points(self):
        start_points = np.empty((self.num_trajs, self.dim))

        for traj_idx, traj in enumerate(self.traj_list):
            start_points[traj_idx] = traj.start_point

        return start_points

    @property
    def end_points(self):
        end_points = np.empty((self.num_trajs, self.dim))

        for traj_idx, traj in enumerate(self.traj_list):
            end_points[traj_idx] = traj.end_point

        return end_points

    """
    Return the total volume of the convex hulls generated by the trajctory pouints at each step.
    Useful for benchmarking strategies against some baseline.
    """
    def total_conv_hull_vol(self):
        conv_hull_data = np.empty((self.traj_len, self.dim))

        for idx in range(self.traj_len):
            conv_hull_data[idx] = ConvexHull(self[idx]).volume

        return np.sum(conv_hull_data)

    """
    Returns 3D matrix labeled by dimensions (step, traj_num, var) such that
    step: the step number of trajectory data
    traj_num: index referring to (traj_num)^{th} trajectory in collection
    var: the variable index
    The combined matrix (traj_pts_mat) essentially contains the trajectory data ordered according to computation step such that
    traj_pts_mat[step] is a 2D matrix such that traj_pts_mat[step][traj_num] contains the (step)^{th} step of the (traj_num)^{th} trajectory in the
    collection.
    """
    def get_combined_step_mat(self):
        traj_pts_mat = np.empty((self.traj_len, self.num_trajs, self.dim))

        for step in range(self.traj_len):
            step_mat = np.empty((self.num_trajs, self.dim))

            for idx, traj in enumerate(self.traj_list):
                step_mat[idx] = traj[step]

            traj_pts_mat[step] = step_mat

        return traj_pts_mat

    def __len__(self):
        return self.num_trajs

    """
    Here, the index refers to the step number of the trajectory.
    Returns a matrix with trajectory data referring to the (index)^{th} step
    """
    def __getitem__(self, index):
        traj_points = np.empty((self.num_trajs, self.dim))

        for idx, traj in enumerate(self.traj_list):
            traj_points[idx] = traj[index]

        return traj_points

    def __iter__(self):
        return iter(self.traj_list)
