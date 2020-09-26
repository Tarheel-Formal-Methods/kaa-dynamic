import random
import sympy as sp
import numpy as np
import multiprocessing as mp

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
        num_traj: number of trajectories to generate.
        time_steps: number of time steps to generate trajs.
@returns list of Traj objects representing num random trajectories.
"""
def generate_traj(bund, num_traj, time_steps):

    model = bund.model
    var = bund.vars

    trajs = [ Traj(model) for _ in range(num_traj) ]

    box_interval = calc_envelop_box(bund)
    initial_points = []
    points_generated = 0

    'Generate points by enveloping a box over the initial polyhedron and picking the points that land within the polyhedron'
    while points_generated < num_traj:
        gen_point = list(map(lambda x: random.uniform(x[0], x[1]), box_interval))
        if check_membership(gen_point, bund):
            initial_points.append(gen_point)
            points_generated += 1

    'Create trajectories with initial points.'
    for point_idx, point in enumerate(initial_points):
        trajs[point_idx].add_point(point)

    df = model.f

    'Parallelize point propagation'
    p = mp.Pool(processes=4)
    prop_trajs = p.starmap(__point_prop_worker, [ (trajs[i], time_steps, df, var) for i in range(num_traj) ])
    p.close()
    p.join()

    return prop_trajs

"""
Worker routine for processes to propagate points for alloted number of time_steps
@params traj: Traj object containing initial point
        time_steps: number of time steps to generate trajectory
        df: transformation dynamics from model
        var: list of variables from model
@returns Traj object containing propagated points.
"""
def __point_prop_worker(traj, time_steps, df, var):
    
    'Propagate the points according to the dynamics for designated number of time steps.'
    prev_point = traj[0]
    
    for _ in range(time_steps):

        var_sub = [ (var, prev_point[var_idx]) for var_idx, var in enumerate(var)]
        next_point = [ f.subs(var_sub) for f in df ]
        traj.add_point(next_point)
        
        prev_point = next_point

    return traj

"""
Calculate the enveloping box over the initial polyhedron
@params model: input model
@returns list of intervals representing edges of box.
"""
def calc_envelop_box(bund):

    A, b = bund.getIntersect()

    box_interval = []
    for i in range(dim):
        y = [0 for _ in range(bund.dim)]
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

"""
Calculates naive supremum bound on the difference between two Flowpipe objects
@params flowpipe1: first Flowpipe object
        flowpipe2: second Flowpipe object
        var_ind: index of variable.
@returns maximum difference calculated along desired projection.
"""
def sup_error_bounds(flowpipe1, flowpipe2, var_ind):
    y1_max, y1_min = flowpipe1.get2DProj(var_ind)
    y2_max, y2_min = flowpipe2.get2DProj(var_ind)

    max_diff = np.absolute(np.subtract(y1_max, y2_max))
    min_diff = np.absolute(np.subtract(y1_min, y2_min))

    return np.amax(np.append(max_diff, min_diff))
