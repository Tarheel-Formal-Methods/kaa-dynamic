from models.phos import Phosphorelay, Phosphorelay_UnitBox
from kaa.templates import MultiStrategy
from kaa.experi_init import *

from kaa.timer import Timer

def test_sapo_Phos():
    num_steps = 100
    model = Phosphorelay()

    experi_input = dict(model=model, #Encompass strat initilizations?
                        strat=None,
                        label="SapoPhos",
                        num_steps=num_steps)

    harosc = ProjectionPlotExperiment(experi_input)
    harosc.execute(0, 1, 2)

def test_sapo_vol_Phos():
    use_supp = True
    use_pregen = False

    num_trajs = 5000
    num_steps = 100

    model = Phosphorelay()
    experi_input = dict(model=model, #Encompass strat initilizations?
                        strat=None,
                        label="SapoVDP",
                        supp_mode = use_supp,
                        pregen_mode = use_pregen,
                        num_trajs=num_trajs,
                        num_steps=num_steps-1,
                        max_steps=num_steps)

    harosc = VolumeExperiment(experi_input)
    harosc.execute(1)

def test_pca_lin_Phos():

    NUM_STEPS = 30
    PHOS_LIN_ITER_STEPS = 1 #Number of steps between each recomputation of LinApp Templates.
    PHOS_PCA_ITER_STEPS = 1 #Number of steps between each recomputation of PCA Templates.
    'PCA Strategy Parameters'
    PHOS_PCA_TRAJ_STEPS = 5 #Number of steps our sample trajectories should run.
    PHOS_PCA_NUM_TRAJ = 200 #Number of sample trajectories we should use for the PCA routine.
    PHOS_LIFE_SPAN = 3

    unit_model = Phosphorelay_UnitBox()
    unit_mod_reach = ReachSet(unit_model)

    #points = [[1,1,1], [1.005, 1,1], [1.01,1.01,1.01], [1.005,1.01,1.01], [1,1.005,1], [1.01,1,1.05]]
    #trajs = [Traj(unit_model , point, NUM_STEPS) for point in points]

    multi_strat = MultiStrategy(LinStrat(unit_model, iter_steps=PHOS_LIN_ITER_STEPS), \
                                DelayedPCAStrat(unit_model, traj_steps=PHOS_PCA_TRAJ_STEPS, num_trajs=PHOS_PCA_NUM_TRAJ, life_span=PHOS_LIFE_SPAN))
                                #PCAStrat(unit_model, traj_steps=PHOS_PCA_TRAJ_STEPS, num_trajs=PHOS_PCA_NUM_TRAJ, iter_steps=PHOS_PCA_ITER_STEPS+PHOS_PCA_DELAY))
                                #PCAStrat(unit_model, traj_steps=PHOS_PCA_TRAJ_STEPS, num_trajs=PHOS_PCA_NUM_TRAJ, iter_steps=PHOS_PCA_ITER_STEPS+2*PHOS_PCA_DELAY))
    mod_lin_flow = unit_mod_reach.computeReachSet(NUM_STEPS, tempstrat=multi_strat)

   # points = [[0,1.97], [0.01, 1.97], [0.01,2], [0,2], [0.005,1.97], [0.005,2], [0,1.97],  [0,1.985], [0.01,1.985]]
    #trajs = [Traj(unit_model, point, NUM_STEPS) for point in points]

    phos_plot = Plot()
    phos_plot.add(mod_lin_flow)

    'Add trajectories'
    for traj in trajs:
        phos_plot.add(traj)

    phos_plot.plot2DPhase(0,1, separate=False, plotvertices=True)

def test_sliding_skewed_plot_Phos():
    use_supp = True
    use_pregen = False

    num_trajs = 5000
    num_steps = 150

    pca_window_size = 10
    lin_window_size = 10

    model = Phosphorelay_UnitBox()

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
    experi.execute(0,1,2)
    Timer.generate_stats()

def test_skewed_sliding_strat_comb_Phos():
    unit_model = Phosphorelay_UnitBox()
    model = Phosphorelay()
    test_skewed_sliding_strat_comb(unit_model, 200, 5000, use_supp=True, use_pregen=False, use_sapo=model)

def test_sliding_pca_Phos():
    model = Phosphorelay_UnitBox()
    test_sliding_pca(model, 20, 200, 5000, use_supp=True, use_pregen=False)

def test_sliding_lin_Phos():
    model = Phosphorelay_UnitBox()
    test_sliding_lin(model, 20, 200, 5000, use_supp=True, use_pregen=False)
