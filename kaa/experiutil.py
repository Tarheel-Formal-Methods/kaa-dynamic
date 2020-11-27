import math
import random as ran
import sympy as sp
import numpy as np
from itertools import product, chain
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
        num_traj: umber of trajectories to generate.
        time_steps: number of time steps to generate trajs.
@returns list of Traj objects representing num random trajectories.
"""
def generate_traj(bund, num_trajs, time_steps):

    model = bund.model

    initial_points = gen_ran_pts_box(bund, num_trajs)

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
Generates random points contained within the tightest enveloping parallelotope of the Chevyshev sphere.
@params bund: Bundle object
        num_trajs: number of trajs to generate
        shrinkfactor: factor to shrink the radius of the sphere. This allows a box with smaller dimensions
@returns list of generated random points.
"""
def gen_ran_pts_box(bund, num_trajs, shrinkfactor=1):
    bund_sys = bund.getIntersect()
    chebycenter = bund_sys.chebyshev_center
    
    center = chebycenter.center
    radius = chebycenter.radius

    box_intervals = [[c - (radius*shrinkfactor), c + (radius*shrinkfactor)] for c in center]

    gen_points = []
    for _ in range(num_trajs):
        gen_points.append([ran.uniform(bound[0], bound[1]) for bound in box_intervals])

    return gen_points

"""
Generates random points through the generators of the initial box and checking membership. Could be extremely time-consuming if the intersection is thin.
@params bund: Bundle object
        num_trajs: number of trajs to generate
@returns points generated.
"""
def gen_ran_pts_ptope(bund, num_trajs):

    bund_sys = bund.getIntersect()
    points_generated = 0
    gen_pts = []

    ptope = bund.ptopes[0] #Initial Box from Bund.
    gen_expr = ptope.getGeneratorRep()

    while points_generated < num_trajs:
        interval_ran_pts = [ (var, ran.uniform(0,1)) for var in bund.vars ]
        ran_pt = [ expr.subs(interval_ran_pts, simultaneous=True) for expr in gen_expr ]

        if bund_sys.check_membership(ran_pt):
            gen_pts.append(ran_pt)
            points_generated += 1

    return gen_pts

"""
Find corner vertices for an initial box along with midpoints between the corners.
@params init_box : intervals of the box given as a list of lists where each member's left,right value
                   are the start,end points respectively for the intervals of the box.
@returns list of border points.
"""
def get_init_box_borders(init_box):

    midpoints = [ start + (end - start) / 2 for start, end in init_box  ]
    border_points = list(product(*init_box))

    for point_idx, point in enumerate(midpoints):
        half_points = [init_inter if point_idx != inter_idx else [point] for inter_idx, init_inter in enumerate(init_box) ]

        border_points += list(product(*half_points))

    return border_points

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
