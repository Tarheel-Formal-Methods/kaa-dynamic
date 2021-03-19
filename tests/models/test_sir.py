from kaa.reach import ReachSet
from kaa.plotutil import Plot
from models.sir import SIR, SIR_UnitBox
from kaa.temp.pca_strat import *
from kaa.temp.lin_app_strat import *
from kaa.templates import MultiStrategy

from kaa.timer import Timer
from kaa.trajectory import Traj
from kaa.experi_init import *

from kaa.bundle import BundleMode

def test_sapo_SIR():
    num_steps = 150
    model = SIR(delta=0.5)

    experi_input = dict(model=model, #Encompass strat initilizations?
                        strat=None,
                        label="SapoSIR",
                        num_steps=num_steps)

    harosc = ProjectionPlotExperiment(experi_input)
    harosc.execute(0,1,2)


    'Generaste the trajectories and add them to the plot.'
    #for traj in trajs:
    #    sir_plot.add(traj)
    sir_plot.add(mod_flow)
    #sir_plot.add(mod_unit_flow)
    sir_plot.plot2DPhase(0,1)

    Timer.generate_stats()

def test_sliding_skewed_plot_SIR():
    use_supp = True
    use_pregen = False

    num_trajs = 5000
    num_steps = 150

    pca_window_size = 10
    lin_window_size = 10

    model = SIR_UnitBox(delta=0.5)

    pca_strat = SlidingPCAStrat(model, lifespan=pca_window_size)
    lin_strat = SlidingLinStrat(model, lifespan=lin_window_size)

    experi_input = dict(model=model,
                        strat=MultiStrategy(pca_strat, lin_strat),
                        label=f"SlidingPCA Step {pca_window_size} and SlidingLin Step {lin_window_size}",
                        supp_mode = use_supp,
                        pregen_mode = use_pregen,
                        num_trajs=num_trajs,
                        num_steps=num_steps-1,
                        max_steps=num_steps)

    if use_supp:
        file_identifier = "(SUPP)"
    elif use_pregen:
        file_identifier = f"(PREGEN: {num_trajs})"
    else:
        file_identifier = "(RAND)"

    experi = ProjectionPlotExperiment(experi_input)
    experi.execute(0, 1, 2)
    Timer.generate_stats()

def test_sir_lin_pca_strat():

    NUM_STEPS = 70
    SIR_PCA_ITER_STEPS = 1 #Number of steps between each recomputation of PCA Templates.
    'PCA Strategy Parameters'
    SIR_PCA_TRAJ_STEPS = 1 #Number of steps our sample trajectories should run.
    SIR_PCA_NUM_TRAJ = 100 #Number of sample trajectories we should use for the PCA routine.
    SIR_LIN_ITER_STEPS = 1
    #
    SIR_PCA_LIFE_SPAN = 3

    sir_pca = SIR_UnitBox(delta=0.5)
    sir_plot = Plot()

    points = [[0.79,0.19,0], [0.79, 0.2,0], [0.8,0.19,0], [0.8,0.2,0], [0.79,0.195,0], [0.8,0.195,0], [0.795,0.19,0],  [0.795,0.2,0]]
    trajs = [Traj(sir_pca, point, NUM_STEPS) for point in points]

    pca_strat = MultiStrategy(LinStrat(sir_pca, iter_steps=SIR_LIN_ITER_STEPS), \
                              DelayedPCAStrat(sir_pca, traj_steps=SIR_PCA_TRAJ_STEPS, num_trajs=SIR_PCA_NUM_TRAJ, life_span=SIR_PCA_LIFE_SPAN))

    sir_pca_reach = ReachSet(sir_pca)
    sir_flow_pca = sir_pca_reach.computeReachSet(NUM_STEPS, tempstrat=pca_strat)
    sir_plot.add(sir_flow_pca, "SIR_LinApp&PCA")

    'Add trajectories'
    for traj in trajs:
        sir_plot.add(traj)

   # sir_plot.plot2DPhase(0,1,separate=False, plotvertices=True)
    sir_plot.plot2DPhase(1, 2, separate=False, plotvertices=True)
    sir_plot.plot2DPhase(0, 2, separate=False, plotvertices=True)

    Timer.generate_stats()

def test_strat_comb_SIR():
    model = SIR_UnitBox(delta=0.5)
    test_strat_comb(model, (1,3,5), 150, -1)

def test_skewed_sliding_strat_comb_SIR():
    model = SIR_UnitBox(delta=0.08)
    test_skewed_sliding_strat_comb(model, 150, 5000, use_supp=True, use_pregen=False)

def test_sliding_strat_comb_SIR():
    model = SIR_UnitBox(delta=0.08)
    test_sliding_strat_comb(model, 150, 5000, use_supp=True, use_pregen=False)

def test_one_one_strat_pca_SIR():
    model = SIR_UnitBox(delta=0.08)
    test_one_one_strat_pca(model, 150)

def test_one_one_strat_lin_SIR():
    model = SIR_UnitBox(delta=0.08)
    test_one_one_strat_lin(model, 150)

def test_sliding_pca_SIR():
    model = SIR_UnitBox(delta=0.5)
    test_sliding_pca(model, 20, 150, 5000, use_supp=True, use_pregen=False)

def test_sliding_lin_SIR():
    model = SIR_UnitBox(delta=0.5)
    test_sliding_lin(model, 20, 150, 5000, use_supp=True, use_pregen=False)

def gen_save_dirs_SIR():
    model = SIR_UnitBox(delta=0.5)
    gen_save_dirs(model, 150)

def find_pca_variation_SIR():
    unit_model = SIR_UnitBox(delta=0.5)
    find_pca_variation(unit_model, 150, max_num_trajs=6000, label="PCADevSIR")

def find_lin_variation_SIR():
    unit_model = SIR_UnitBox(delta=0.5)
    find_lin_variation(unit_model, 150, max_num_trajs=6000, label="LinDevSIR")
