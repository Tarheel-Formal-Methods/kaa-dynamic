import random
import sympy as sp
import numpy as np
import multiprocessing as mp

from kaa.trajectory import Traj
from kaa.settings import KaaSettings


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
    bund_sys = bund.getIntersect()
 
    points_generated = 0

    'Choose the first bundle to generate random points over.'
    ptope = bund.ptopes[0]
    initial_points = []

    while points_generated < num_traj:
        ran_pt = ptope.gen_random_pt()
        if bund_sys.check_membership(ran_pt):
            initial_points.append(ran_pt)
            points_generated += 1

    trajs = [ Traj(model, point, steps=time_steps) for point in initial_points ]

    """
    if KaaSettings.use_parallel:
        'Parallelize point propagation'
        p = mp.Pool(processes=4)
        prop_trajs = p.starmap(point_prop, [ (model, point, time_steps) for point in initial_points ])
        p.close()
        p.join()
    """

    return trajs

"""
Calculate the enveloping box over the initial polyhedron
@params model: input model
@returns list of intervals representing edges of box.
"""
def calc_envelop_box(bund):

    bund_sys = bund.getIntersect()
    box_interval = []

    for i in range(bund.dim):

        y = [0 for _ in range(bund.dim)]
        y[i] = 1
        
        maxCood = bund_sys.max_opt(y).fun
        minCood = bund_sys.min_opt(y).fun
        box_interval.append([minCood, maxCood])

    return box_interval

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
