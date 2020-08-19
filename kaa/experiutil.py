import random
import sympy as sp

from kaa.trajectory import Traj


"""
Generate random trajectories from initial set of model.
@params model: Model
        num: number of trajectories to generate.
        time_steps: number of time steps to generate trajs.
@returns list of Traj objects representing num random trajectories.
"""

def generate_traj(model, num, time_steps):

    assert model.initial_set is not None, "Initial set must be specified in Kaa.Model"

    vars = model.vars

    trajs = [ Traj(model) for _ in range(num) ]
    initial_points = [ list(map(lambda x: random.uniform(x[0], x[1]), model.initial_set)) for _ in range(num) ]

    'Create trajectories with initial points.'
    for point_idx, point in enumerate(initial_points):
        trajs[point_idx].add_traj_point(point)

    df = model.f

    prev_points = initial_points
    for _ in range(time_steps):

        next_points_list = []
        for point_idx, point in enumerate(prev_points):
            var_sub = [ (var, point[var_idx]) for var_idx, var in enumerate(vars)]
            next_point = [ f.subs(var_sub) for f in df ]
            trajs[point_idx].add_traj_point(next_point)

            next_points_list.append(next_point)

        prev_points = next_points_list

    return trajs
