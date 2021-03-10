import math
import random as rand
import sympy as sp
import numpy as np
from itertools import product, chain

from kaa.temp.pca_strat import *
from kaa.temp.lin_app_strat import *
from kaa.templates import MultiStrategy
from kaa.experiment import *
from kaa.settings import KaaSettings
from kaa.log import Output

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

"""
def test_comp_ani(model, x, y, num_steps)

    NUM_STEPS = 70
    #model = VanDerPol(delta=0.08)
    unit_model = VanDerPol_UnitBox(delta=0.08)

    VDP_PCA_ITER_STEPS = 1 #Number of steps between each recomputation of PCA Templates.
    'PCA Strategy Parameters'
    VDP_PCA_TRAJ_STEPS = 5 #Number of steps our sample trajectories should run.
    VDP_PCA_NUM_TRAJ = 100 #Number of sample trajectories we should use for the PCA routine.
    VDP_PCA_DELAY = 5

    pca_dirs = GeneratedPCADirs(unit_model, VDP_PCA_NUM_TRAJ, NUM_STEPS)
    pca_strat = PCAStrat(unit_model, traj_steps=VDP_PCA_TRAJ_STEPS, num_trajs=VDP_PCA_NUM_TRAJ, iter_steps=VDP_PCA_ITER_STEPS, pca_dirs=pca_dirs)
    experi_input = dict(model=unit_model,
                   strat=pca_strat,
                   label="",
                   num_steps=NUM_STEPS)

    vdp_pca = Animation(experi_input)
    vdp_pca.execute()
    vdp_pca.animate(0,1, pca_strat)
"""

def test_sliding_strat_comb(model, num_steps, num_trajs, num_trials=10, use_supp=True, use_pregen=False):
    if use_supp:
        num_trials = 1

    inputs = []
    for window_size in range(10,0,-2):
        pca_strat = SlidingPCAStrat(model, lifespan=window_size)
        lin_strat = SlidingLinStrat(model, lifespan=window_size)

        experi_input = dict(model=model,
                            strat=MultiStrategy(pca_strat, lin_strat),
                            label=f"SlidingPCA Step {window_size} and SlidingLin Step {window_size}",
                            supp_mode = use_supp,
                            pregen_mode = use_pregen,
                            num_trajs=num_trajs,
                            num_steps=num_steps-1,
                            max_steps=num_steps)

        inputs.append(experi_input)

    if use_supp:
        file_identifier = "(SUPP)"
    elif use_pregen:
        file_identifier = "(PREGEN: {num_trajs})"
    else:
        file_identifier = "(RAND)"

    experi = VolumeExperiment(*inputs, label=f"SlidingCombination {model.name} {file_identifier}")
    experi.execute(num_trials)

def test_strat_comb(model, step_tup, num_steps, num_trajs, num_trials=10, use_supp=False, use_pregen=True):
    if use_supp:
        num_trials = 1

    pca_iter_steps = step_tup
    lin_iter_steps = step_tup

    inputs = []
    for pca_step, lin_step in product(pca_iter_steps, lin_iter_steps): #model tossed around too many times.
        pca_strat = PCAStrat(model, iter_steps=pca_step)
        lin_strat = LinStrat(model, iter_steps=lin_step)
        experi_input = dict(model=model,
                            strat=MultiStrategy(pca_strat, lin_strat),
                            label=f"PCA Step {pca_step} and Lin Step {lin_step}",
                            supp_mode = use_supp,
                            pregen_mode = use_pregen,
                            num_trajs=num_trajs,
                            num_steps=num_steps-1,
                            max_steps=num_steps)

        inputs.append(experi_input)

    if use_supp:
        file_identifier = "(SUPP)"
    elif use_pregen:
        file_identifier = "(PREGEN: {num_trajs})"
    else:
        file_identifier = "(RAND)"

    experi = VolumeExperiment(*inputs, label=f"Combination with PCA and Lin Strats {model.name} {file_identifier}")
    experi.execute(num_trials)

def test_one_one_strat_pca(model, num_steps, num_trajs, num_trials=10, use_supp=False, use_pregen=True):
    if use_supp:
        num_trials = 1

    inputs = []
    for pca_step in range(2,6): #model tossed around too many times.
        pca_strat1 = PCAStrat(model, iter_steps=1)
        lin_strat1 = LinStrat(model, iter_steps=1)
        experi_strat1 = PCAStrat(model, iter_steps=pca_step)
        experi_input1 = dict(model=model,
                             strat=MultiStrategy(pca_strat1, lin_strat1, experi_strat1),
                             label=f"PCA with Box and 1-1 Temp for {pca_step} steps",
                             supp_mode = use_supp,
                             pregen_mode = use_pregen,
                             num_trajs=num_trajs,
                             num_steps=num_steps-1,
                             max_steps=num_steps)

        inputs.append(experi_input1)

    if use_supp:
        file_identifier = "(SUPP)"
    elif use_pregen:
        file_identifier = "(PREGEN: {num_trajs})"
    else:
        file_identifier = "(RAND)"

    experi = VolumeExperiment(*inputs, label=f"1-1 Strat Base PCA Trials {model.name} {file_identifier}")
    experi.execute(num_trials)

def test_one_one_strat_lin(model, num_steps, num_trajs, num_trials=10, use_supp=False, use_pregen=True):
    if use_supp:
        num_trials = 1

    inputs = []
    for lin_step in range(2,max_step+1): #model tossed around too many times.
        pca_strat1 = PCAStrat(model, iter_steps=1)
        lin_strat1 = LinStrat(model, iter_steps=1)
        experi_strat1 = LinStrat(model, iter_steps=lin_step)
        experi_input1 = dict(model=model,
                            strat=MultiStrategy(pca_strat1, lin_strat1, experi_strat1),
                            label=f"PCA with Box and 1-1 Temp for {lin_step} steps",
                             supp_mode = use_supp,
                             pregen_mode = use_pregen,
                             num_trajs=num_trajs,
                             num_steps=num_steps-1,
                             max_steps=num_steps)

        inputs.append(experi_input1)

    if use_supp:
        file_identifier = "(SUPP)"
    elif use_pregen:
        file_identifier = "(PREGEN: {num_trajs})"
    else:
        file_identifier = "(RAND)"

    experi = VolumeExperiment(*inputs, label=f"1-1 Strat Base LinApp Trials {model.name} {file_identifier}")
    experi.execute(num_trials)

def test_sliding_pca(model, max_life, num_steps, num_trajs, life_incre=5, num_trials=10 , use_supp=False, use_pregen=True):
    if use_supp:
        num_trials = 1

    NUM_STEPS = num_steps
    NUM_TRAJ = num_trajs #Number of sample trajectories we should use for the PCA routine.
    LIFE_MAX = max_life
    LIFE_INCREMENT = life_incre

    inputs = []
    for lifespan in range(LIFE_MAX, 0, -LIFE_INCREMENT): #model tossed around too many times.
        experi_strat = SlidingPCAStrat(model, lifespan=lifespan)
        experi_input = dict(model=model,
                            strat=experi_strat,
                            label=f"SlidingPCA with Window Size:{lifespan}",
                            supp_mode = use_supp,
                            pregen_mode = use_pregen,
                            num_trajs=NUM_TRAJ,
                            num_steps=NUM_STEPS-1,
                            max_steps=NUM_STEPS)
        inputs.append(experi_input)

    for lifespan in range(LIFE_INCREMENT, 0, -1): #model tossed around too many times.
        experi_strat = SlidingPCAStrat(model, lifespan=lifespan)
        experi_input = dict(model=model,
                            strat=experi_strat,
                            label=f"SlidingPCA with Window Size:{lifespan}",
                            supp_mode = use_supp,
                            pregen_mode = use_pregen,
                            num_trajs=NUM_TRAJ,
                            num_steps=NUM_STEPS-1,
                            max_steps=NUM_STEPS)

        inputs.append(experi_input)

    if use_supp:
        file_identifier = "(SUPP)"
    elif use_pregen:
        file_identifier = "(PREGEN: {num_trajs})"
    else:
        file_identifier = "(RAND)"

    experi = VolumeExperiment(*inputs, label=f"SlidingPCA{model.name} {file_identifier}")
    experi.execute(num_trials)

def test_sliding_lin(model, max_life, num_steps, num_trajs, life_incre=5, num_trials=10, use_supp=False, use_pregen=True):
    if use_supp:
        num_trials = 1

    NUM_STEPS = num_steps
    LIFE_MAX = max_life
    NUM_TRAJS = num_trajs
    LIFE_INCREMENT = life_incre

    inputs = []
    for lifespan in range(LIFE_MAX, 0, -LIFE_INCREMENT): #model tossed around too many times.
        experi_strat = SlidingLinStrat(model, lifespan=lifespan)
        experi_input = dict(model=model,
                            strat=experi_strat,
                            label=f"SlidingLin with Window Size:{lifespan}",
                            supp_mode = use_supp,
                            pregen_mode = use_pregen,
                            num_trajs=NUM_TRAJS,
                            num_steps=NUM_STEPS-1,
                            max_steps=NUM_STEPS)

        inputs.append(experi_input)

    for lifespan in range(LIFE_INCREMENT, 0, -1): #model tossed around too many times.
        experi_strat = SlidingLinStrat(model, lifespan=lifespan)
        experi_input = dict(model=model,
                            strat=experi_strat,
                            label=f"SlidingLin with Window Size:{lifespan}",
                            supp_mode = use_supp,
                            pregen_mode = use_pregen,
                            num_trajs=NUM_TRAJS,
                            num_steps=NUM_STEPS-1,
                            max_steps=NUM_STEPS)

        inputs.append(experi_input)

    if use_supp:
        file_identifier = "(SUPP)"
    elif use_pregen:
        file_identifier = "(PREGEN: {num_trajs})"
    else:
        file_identifier = "(RAND)"

    experi = VolumeExperiment(*inputs, label=f"SlidingLin{model.name} {file_identifier}")
    experi.execute(num_trials)

def test_comb_stdev_reduction(model, num_steps, num_trials=10):
    if use_supp:
        num_trials = 1

    NUM_STEPS = num_steps

    inputs = []
    for num_trajs in range(1000,4000,1000): #model tossed around too many times.
        pca_strat = PCAStrat(model, iter_steps=1)
        lin_strat = LinStrat(model, iter_steps=1)
        experi_input = dict(model=model,
                            strat=MultiStrategy(pca_strat, lin_strat),
                            label=f"PCA Step 1 and Lin Step 1 with NUM_TRAJ:{num_trajs}",
                            pregen_mode = True,
                            num_trajs=num_trajs,
                            num_steps=NUM_STEPS,
                            max_steps=NUM_STEPS)

        inputs.append(experi_input)

    experi = VolumeExperiment(*inputs, label="Combination with PCA and Lin Strats")
    experi.execute(num_trials)

def gen_save_dirs(model, num_steps, max_num_trajs=8000, num_trials=10):
    for num_trajs in range(1000, max_num_trajs+1000, 1000):
        generated_dirs = []
        for trial_num in range(num_trials):
            Output.prominent(f"GENERATED DIRECTIONS FOR TRIAL {trial_num} WITH {num_trajs} TRAJS FOR {num_steps} STEPS")
            gen_pca_dirs = GeneratedPCADirs(model, num_steps, num_trajs)
            gen_lin_dirs = GeneratedLinDirs(model, num_steps, num_trajs)

            gen_dirs_tuple = GenDirsTuple(gen_pca_dirs, gen_lin_dirs)
            generated_dirs.append(gen_dirs_tuple)

            update_seed()

        reset_seed()
        DirSaveLoader.save_dirs(model, num_steps, num_trajs, KaaSettings.RandSeed, generated_dirs)

def find_pca_variation(model, num_steps, num_trials=10, max_num_trajs=8000, label=""):
    inputs = []
    for num_trajs in range(1000, max_num_trajs+1000, 1000):
        experi_input = dict(model=model,
                            strat=None,
                            label= f"PCA from {num_trajs} TRAJS",
                            num_trajs=num_trajs,
                            pregen_mode = True,
                            num_steps=num_steps,
                            max_steps=num_steps)

        inputs.append(experi_input)

    experi = DeviationExperiment(*inputs, "PCADev", label=label)
    experi.execute(num_trials)

def find_lin_variation(model, num_steps, num_trials=1, max_num_trajs=8000, label=""):
    inputs = []
    for num_trajs in range(1000, max_num_trajs+1000, 1000):
        experi_input = dict(model=model,
                            strat=None,
                            label= f"LinApp from {num_trajs} TRAJS",
                            pregen_mode = True,
                            num_trajs=num_trajs,
                            num_steps=num_steps,
                            max_steps=num_steps)
        inputs.append(experi_input)

    experi = DeviationExperiment(*inputs, "LinDev", label=label)
    experi.execute(num_trials)
