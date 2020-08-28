import random
import sympy as sp
import numpy as np

from kaa.trajectory import Traj
from kaa.lputil import minLinProg, maxLinProg

"""
Generate random trajectories from initial set (initial bundle) of model.
@params model: Model
        num: number of trajectories to generate.
        time_steps: number of time steps to generate trajs.
@returns list of Traj objects representing num random trajectories.
"""
def generate_init_traj(model, num, time_steps):

    bund = model.bund

    return generate_traj(bund, num, time_steps)

"""
Generate random trajectories from polytope defined by parallelotope bundle.
@params model: Model
        num: number of trajectories to generate.
        time_steps: number of time steps to generate trajs.
@returns list of Traj objects representing num random trajectories.
"""
def generate_traj(bund, num, time_steps):

    model = bund.model
    vars = bund.vars

    trajs = [ Traj(model) for _ in range(num) ]

    box_interval = calc_envelop_box(bund)
    initial_points = []
    points_generated = 0

    'Generate points by enveloping a box over the initial polyhedron and picking the points that land within the polyhedron'
    while points_generated < num:
        gen_point = list(map(lambda x: random.uniform(x[0], x[1]), box_interval))
        if check_membership(gen_point, bund):
            initial_points.append(gen_point)
            points_generated += 1

    'Create trajectories with initial points.'
    for point_idx, point in enumerate(initial_points):
        trajs[point_idx].add_traj_point(point)

    df = model.f

    'Propagate the points according to the dynamics for designated number of time steps.'
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

"""
Calculate the enveloping box over the initial polyhedron
@params model: input model
@returns list of intervals representing edges of box.
"""
def calc_envelop_box(bund):

    A, b = bund.getIntersect()
    dim = len(bund.vars)

    box_interval = []
    for i in range(dim):
        y = [0 for _ in range(dim)]
        y[i] = 1
        
        maxCood = maxLinProg(y, A, b).fun
        minCood = minLinProg(y, A, b).fun
        box_interval.append([minCood, maxCood])

    return box_interval
"""
Checks if point is indeed contained in initial polyhedron
@params point: point to test
        model: input model
@returns boolean value indictating membership.
"""
def check_membership(point, bund):

    A, b = bund.getIntersect()

    for row_idx, row in enumerate(A):
        if np.dot(row, point) > b[row_idx]:
            return False
    return True
