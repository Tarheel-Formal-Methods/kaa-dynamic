from kaa.reach import ReachSet
from kaa.plotutil import Plot
from models.quadcopter import Quadcopter, Quadcopter_UnitBox
from kaa.experi_init import *

from kaa.timer import Timer

def test_sapo_Quad():

    num_steps = 300
    model = Quadcopter()

    experi_input = dict(model=model, #Encompass strat initilizations?
                        strat=None,
                        label="SapoQuad",
                        num_steps=num_steps)

    harosc = ProjectionPlotExperiment(experi_input)
    harosc.execute(2,6,17)

def test_sliding_skewed_plot_Quad():
    use_supp = True
    use_pregen = False

    num_trajs = 5000
    num_steps = 150

    pca_window_size = 10
    lin_window_size = 10

    model = Quadcopter_UnitBox()

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
    experi.execute(2,6,17)
    Timer.generate_stats()

"""
def test_pca_Quad():
    NUM_STEPS = 5
    QUAD_PCA_ITER_STEPS = 1 #Number of steps between each recomputation of PCA Templates.
    'PCA Strategy Parameters'
    QUAD_PCA_TRAJ_STEPS = 1 #Number of steps our sample trajectories should run.
    QUAD_PCA_NUM_TRAJ = 200 #Number of sample trajectories we should use for the PCA routine.

    quad_pca = Quadcopter_UnitBox()
    quad_plot = Plot()

    #trajs = [Traj(rossler_pca, point, NUM_STEPS) for point in points]

    pca_strat = PCAStrat(quad_pca, traj_steps=QUAD_PCA_TRAJ_STEPS, num_trajs=QUAD_PCA_NUM_TRAJ, iter_steps=QUAD_PCA_ITER_STEPS)
    quad_pca_reach = ReachSet(quad_pca)
    quad_flow_pca = quad_pca_reach.computeReachSet(NUM_STEPS, tempstrat=pca_strat)

    'Add trajectories'
    #for traj in trajs:
    #    rossler_plot.add(traj)

    quad_plot.add(quad_flow_pca, "Quad_PCA")
    quad_plot.plot2DPhase(2,5, separate=True, plotvertices=True)
"""

def test_skewed_sliding_strat_comb_Quad():
    model = Quadcopter_UnitBox()
    test_skewed_sliding_strat_comb(model, 150, 5000, use_supp=True, use_pregen=False)

def test_sliding_pca_Quad():
    model = Quadcopter_UnitBox()
    test_sliding_pca(model, 20, 150, 5000, use_supp=True, use_pregen=False)

def test_sliding_lin_Quad():
    model = Quadcopter_UnitBox()
    test_sliding_lin(model, 20, 150, 5000, use_supp=True, use_pregen=False)
