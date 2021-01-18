import math
import random as rand
import sympy as sp
import numpy as np
from itertools import product, chain

from kaa.temp.pca_strat import *
from kaa.temp.lin_app_strat import *
from kaa.templates import MultiStrategy
from kaa.experiment import Experiment, ExperimentBatch, exec_plot_vol_results
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
    return bund.generate_traj(num, time_steps)

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
        interval_ran_pts = [(var, rand.uniform(0,1)) for var in bund.vars]
        ran_pt = [expr.subs(interval_ran_pts, simultaneous=True) for expr in gen_expr]

        if bund_sys.check_membership(ran_pt):
            gen_pts.append(ran_pt)
            points_generated += 1

    return gen_pts

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

def test_strat_comb(model, step_tup, num_steps, num_trajs, num_trials=10, filename="STRATCOMB"):
    NUM_STEPS = num_steps
    NUM_TRAJ = num_trajs #Number of sample trajectories we should use for the PCA routine.
    MAX_STEP = max(step_tup)

    pca_iter_steps = step_tup
    lin_iter_steps = step_tup

    inputs = []
    for pca_step, lin_step in product(pca_iter_steps, lin_iter_steps): #model tossed around too many times.
        pca_strat = PCAStrat(model, iter_steps=pca_step)
        lin_strat = LinStrat(model, iter_steps=lin_step)
        experi_input = dict(model=model,
                            strat=MultiStrategy(pca_strat, lin_strat),
                            label=f"PCA Step {pca_step} and Lin Step {lin_step}",
                            num_trajs=NUM_TRAJ,
                            num_steps=NUM_STEPS)
        inputs.append(experi_input)

    experi = Experiment(*inputs, label="Combination with PCA and Lin Strats", num_trials=num_trials)
    experi.execute()

def test_one_one_strat_pca(model, max_step, num_steps, num_trials=10, filename="ONEONEPCA"):
    NUM_STEPS = num_steps
    NUM_TRAJ = 300 #Number of sample trajectories we should use for the PCA routine.
    LIN_ITER_STEPS = 1 #Number of steps between each recomputation of LinApp Templates.
    PCA_ITER_STEPS = 1 #Number of steps between each recomputation of PCA Templates.
    MAX_STEP = max_step

    inputs = []
    for pca_step in range(2,6): #model tossed around too many times.
        pca_strat1 = PCAStrat(model, iter_steps=1)
        lin_strat1 = LinStrat(model, iter_steps=1)
        experi_strat1 = PCAStrat(model, iter_steps=pca_step)
        experi_input1 = dict(model=model,
                             strat=MultiStrategy(pca_strat1, lin_strat1, experi_strat1),
                             label=f"PCA with Box and 1-1 Temp for {pca_step} steps",
                             num_trajs=NUM_TRAJ,
                             num_steps=NUM_STEPS)
        inputs.append(experi_input1)

    experi = Experiment(*inputs, label="1-1 Strat Base PCA Trials", num_trials=num_trials)
    experi.execute()

def test_one_one_strat_lin(model, max_step, num_steps, num_trials=10, filename="ONEONELIN"):
    NUM_STEPS = num_steps
    NUM_TRAJ = 1000 #Number of sample trajectories we should use for the PCA routine.
    MAX_STEP = max_step

    inputs = []
    for lin_step in range(2,max_step+1): #model tossed around too many times.
        pca_strat1 = PCAStrat(model, iter_steps=1)
        lin_strat1 = LinStrat(model, iter_steps=1)
        experi_strat1 = LinStrat(model, iter_steps=lin_step)
        experi_input1 = dict(model=model,
                            strat=MultiStrategy(pca_strat1, lin_strat1, experi_strat1),
                            label=f"PCA with Box and 1-1 Temp for {lin_step} steps",
                            num_trajs=NUM_TRAJ,
                            num_steps=NUM_STEPS)
        inputs.append(experi_input1)

    experi = Experiment(*inputs, label="1-1 Strat Base LinApp Trials", num_trials=num_trials)
    experi.execute()

def test_sliding_pca(model, max_life, num_steps, life_incre=5, num_trials=10, filename="SLIDINGPCA"):
    NUM_STEPS = num_steps
    NUM_TRAJ = 1000 #Number of sample trajectories we should use for the PCA routine.
    LIFE_MAX = max_life
    LIFE_INCREMENT = life_incre


    inputs = []
    """
    for lifespan in range(LIFE_MAX, 0, -LIFE_INCREMENT): #model tossed around too many times.
        experi_strat = SlidingPCAStrat(model, lifespan=lifespan)
        experi_input = dict(model=model,
                            strat=experi_strat,
                            label=f"SlidingPCA with Window Size:{lifespan}",
                            num_trajs=NUM_TRAJ,
                            num_steps=NUM_STEPS)
        inputs.append(experi_input)
    """
    for lifespan in range(1, 0, -1): #model tossed around too many times.
        experi_strat = SlidingPCAStrat(model, lifespan=lifespan)
        experi_input = dict(model=model,
                            strat=experi_strat,
                            label=f"SlidingPCA with Window Size:{lifespan}",
                            num_trajs=NUM_TRAJ,
                            num_steps=NUM_STEPS)
        inputs.append(experi_input)


    experi = Experiment(*inputs, label=f"UpperSlidingPCA with NUM_TRAJ:{NUM_TRAJ}", num_trials=num_trials)
    experi.execute()

def test_sliding_lin(model, max_life, num_steps, num_trials=10, life_incre=5, filename="SLIDINGLIN"):
    NUM_STEPS = num_steps
    NUM_TRAJ = 1000 #Number of sample trajectories we should use for the PCA routine.
    LIFE_MAX = max_life
    LIFE_INCREMENT = life_incre

    inputs = []
    for lifespan in range(LIFE_MAX, 0, -LIFE_INCREMENT): #model tossed around too many times.
        experi_strat = SlidingLinStrat(model, lifespan=lifespan)
        experi_input = dict(model=model,
                            strat=experi_strat,
                            label=f"SlidingLin with Window Size:{lifespan}",
                            num_trajs=NUM_TRAJ,
                            num_steps=NUM_STEPS)
        inputs.append(experi_input)

    for lifespan in range(LIFE_INCREMENT - 1, 0, -1): #model tossed around too many times.
        experi_strat = SlidingLinStrat(model, lifespan=lifespan)
        experi_input = dict(model=model,
                            strat=experi_strat,
                            label=f"SlidingLin with Window Size:{lifespan}",
                            num_trajs=NUM_TRAJ,
                            num_steps=NUM_STEPS)
        inputs.append(experi_input)

    experi = Experiment(*inputs, label=f"SlidingPCA with NUM_TRAJ:{NUM_TRAJ}", num_trials=num_trials)
    experi.execute()


def test_comb_stdev_reduction(model, num_steps, num_trials=10, filename="STRATCOMBSTDEV"):
    NUM_STEPS = num_steps
    MAX_STEP = 1
    PCA_ITER_STEPS = 1
    LIN_ITER_STEPS = 1

    inputs = []
    for num_trajs in range(1000,4000,1000): #model tossed around too many times.
        pca_strat = PCAStrat(model, iter_steps=PCA_ITER_STEPS)
        lin_strat = LinStrat(model, iter_steps=LIN_ITER_STEPS)
        experi_input = dict(model=model,
                            strat=MultiStrategy(pca_strat, lin_strat),
                            label=f"PCA Step 1 and Lin Step 1 with NUM_TRAJ:{num_trajs}",
                            num_trajs=num_trajs,
                            num_steps=NUM_STEPS)

        inputs.append(experi_input)

    experi = Experiment(*inputs, label="Combination with PCA and Lin Strats", num_trials=num_trials)
    experi.execute()
