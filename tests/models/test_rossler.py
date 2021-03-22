from kaa.reach import ReachSet
from kaa.plotutil import Plot
from models.rossler import Rossler, Rossler_UnitBox
from kaa.trajectory import Traj
from kaa.experi_init import *
from kaa.timer import Timer

def test_sapo_Rossler():

    model = Rossler()
    num_steps = 150

    experi_input = dict(model=model, #Encompass strat initilizations?
                        strat=None,
                        label="SapoRossler",
                        num_steps=num_steps)

    harosc = ProjectionPlotExperiment(experi_input)
    harosc.execute(0,1,2)

def test_rossler_phase():

    model = Rossler()
    mod_reach = ReachSet(model)
    mod_flow = mod_reach.computeReachSet(200)

    rossler_plot = Plot()
    rossler_plot.add(mod_flow)
    rossler_plot.plot2DPhase(0,1)

    Timer.generate_stats()


def test_sliding_skewed_plot_Rossler():
    use_supp = True
    use_pregen = False

    num_trajs = 5000
    num_steps = 150

    pca_window_size = 10
    lin_window_size = 10

    model = Rossler_UnitBox(delta=0.5)

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


def test_pca_lin_Rossler():
    NUM_STEPS = 5
    ROSS_PCA_ITER_STEPS = 1 #Number of steps between each recomputation of PCA Templates.
    'PCA Strategy Parameters'
    ROSS_PCA_TRAJ_STEPS = 1 #Number of steps our sample trajectories should run.
    ROSS_PCA_NUM_TRAJ = 200 #Number of sample trajectories we should use for the PCA routine.

    rossler_pca = Rossler_UnitBox(delta=0.5)
    rossler_plot = Plot()

    points = [[0.05,4.95,0.05], [0.1,4.95,0.05], [0.05,5,0.05], [0.1,5,0.05], [0.05,4.95,0.05], [0.05,4.95,0.1], [0.1,4.95,0.1], [0.1,5,0.1]]
    trajs = [Traj(rossler_pca, point, NUM_STEPS) for point in points]

    pca_strat = PCALinStrat(rossler_pca, traj_steps=ROSS_PCA_TRAJ_STEPS, num_trajs=ROSS_PCA_NUM_TRAJ, iter_steps=ROSS_PCA_ITER_STEPS)

    ross_pca_reach = ReachSet(rossler_pca)
    ross_flow_pca = ross_pca_reach.computeReachSet(NUM_STEPS, tempstrat=pca_strat)
    rossler_plot.add(ross_flow_pca, "SIR_LinApp&PCA")

    'Add trajectories'
    for traj in trajs:
        rossler_plot.add(traj)

    rossler_plot.plot2DPhase(0,1,separate=True, plotvertices=True)
    Timer.generate_stats()

def test_strat_comb_Rossler():
    model = Rossler_UnitBox()
    test_strat_comb(model, (1,3,5), 150, )

def test_sliding_strat_comb_Rossler():
    model = Rossler_UnitBox()
    test_sliding_strat_comb(model, 150, 4000, use_supp=True, use_pregen=False)

def test_skewed_sliding_strat_comb_Rossler():
    model = Rossler_UnitBox()
    test_skewed_sliding_strat_comb(model, 150, 4000, use_supp=True, use_pregen=False)

def test_sliding_pca_Rossler():
    model = Rossler_UnitBox()
    test_sliding_pca(model, 20, 150, -1, use_supp=True, use_pregen=False)

def test_sliding_lin_Rossler():
    model = Rossler_UnitBox()
    test_sliding_lin(model, 20, 150, -1, use_supp=True, use_pregen=False)

def gen_save_dirs_Rossler():
    model = Rossler_UnitBox()
    gen_save_dirs(model, 150)
