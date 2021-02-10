import numpy as np
from scipy.spatial import ConvexHull

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
            var_sub = [(var, prev_point[var_idx]) for var_idx, var in enumerate(self.vars)]
            next_point = [f.subs(var_sub) for f in df]
            self.add_point(next_point)

            prev_point = next_point

    """
    Returns the projection of the trajectory onto an variable axis
    @params var: var to project onto.
    @returns projection onto axis determined by var.
    """
    def get_proj(self, var):
        return self.traj_set[var]

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
        return [(self.traj_set[var])[index] for var in self.vars]

    def __len__(self):
        return self.num_points

class TrajCollection:

    def __init__(self, model, traj_list):
        assert isinstance(traj_list, list), "input must be list of Traj objects"
        self.traj_list = traj_list
        self.num_trajs = len(self.traj_list)
        self.model = model

    @property
    def start_points(self):
        start_points_list = [traj.start_point for traj in self.traj_list]
        return np.asarray(start_points_list)

    @property
    def end_points(self):
        end_points_list = [traj.end_point for traj in self.traj_list]
        return np.asarray(end_points_list)

    @property
    def max_traj_len(self):
        return max([len(traj) for traj in self.traj_list])

    @property
    def conv_hull_vol(self):
        vol_data = []
        for i in range(self.max_traj_len):
            vol_data.append(ConvexHull([traj[i] for traj in self.traj_list]).volume)

        return sum(vol_data)

    def get_mat(self):
        traj_pts_mat = []

        for step in range(self.max_traj_len):
            step_mat = np.empty((self.num_trajs, self.model.dim))

            for idx, traj in enumerate(self.traj_list):
                step_mat[idx] = traj[step]

            traj_pts_mat.append(step_mat)

        return traj_pts_mat

    def __len__(self):
        return self.num_trajs

    def __getitem__(self, index):
        traj_points_list = [traj[index] for traj in self.traj_list]
        return np.asarray(traj_points_list)

    def __iter__(self):
        return iter(self.traj_list)
