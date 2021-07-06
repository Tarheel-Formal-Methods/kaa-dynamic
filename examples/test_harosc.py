from models.harosc import HarOsc
from kaa.temp.pca_strat import *
from kaa.temp.lin_app_strat import *
from kaa.templates import MultiStrategy
from kaa.experi_lib import *

from settings import PlotSettings
from itertools import product

PlotSettings.save_fig = False

def test_HarOsc():
    num_steps = 10
    model = HarOsc()

    experi_input = dict(model=model, #Encompass strat initilizations?
                            strat=None,
                            label="HarOsc",
                            num_steps=num_steps)

    harosc = ProjectionPlotExperiment(experi_input)
    harosc.execute(0,1)

def test_sliding_lin_HarOsc():
    model = HarOsc()

    lin_strat = LinStrat(model, iter_steps=1)
    experi_input = dict(model=model, #Encompass strat initilizations?
                            strat=lin_strat,
                            label="HarOsc",
                            num_steps=10,
                            max_steps=5,
                            num_trajs=10,
                            supp_mode=True,
                            pregen_mode=False)

    harosc = ProjectionPlotExperiment(experi_input)
    harosc.execute(0,1)

def test_lin_HarOsc():
    model = HarOsc()

    lin_strat = LinStrat(model, iter_steps=1)
    experi_input = dict(model=model, #Encompass strat initilizations?
                            strat=lin_strat,
                            label="HarOsc",
                            num_steps=5,
                            max_steps=5,
                            num_trajs=10,
                            supp_mode=False,
                            pregen_mode=False)

    harosc = PhasePlotExperiment(experi_input)
    harosc.execute(0,1, separate=True)

def test_pca_lin_HarOsc():

    NUM_STEPS = 4
    model = HarOsc()
    #trajs = generate_traj(model, 10, 200)
    mod_reach = ReachSet(model)
    #mod_flow = mod_reach.computeReachSet()

    sir_plot = Plot()

    SIR_LIN_ITER_STEPS = 1 #Number of steps between each recomputation of PCA Templates.
    SIR_PCA_ITER_STEPS = 1 #Number of steps between each recomputation of PCA Templates.
    'PCA Strategy Parameters'
    SIR_PCA_TRAJ_STEPS = 2 #Number of steps our sample trajectories should run.
    SIR_PCA_NUM_TRAJ = 100 #Number of sample trajectories we should use for the PCA routine.

    tandem_strat = MultiStrategy(LinStrat(model, iter_steps=SIR_LIN_ITER_STEPS), PCAStrat(model, traj_steps=SIR_PCA_TRAJ_STEPS, num_trajs=SIR_PCA_NUM_TRAJ, iter_steps=SIR_PCA_ITER_STEPS))
    mod_pca_flow = mod_reach.computeReachSet(NUM_STEPS, tempstrat=tandem_strat, transmode=BundleMode.AFO)
    #trajs = generate_traj(model, 10, 200)

    'Generaste the trajectories and add them to the plot.'
    sir_plot.add(mod_pca_flow, "HarOsc PCA")
    sir_plot.plot2DPhase(0,1)

def test_sliding_pca_HarOsc():
    unit_model = HarOsc()

    NUM_STEPS = 5
    VDP_PCA_NUM_TRAJ = 300 #Number of sample trajectories we should use for the PCA routine.
    SIR_LIN_ITER_STEPS = 1 #Number of steps between each recomputation of PCA Templates.
    pca_dirs = GeneratedPCADirs(unit_model, VDP_PCA_NUM_TRAJ, NUM_STEPS+1) #Way to deduce lengeth beforehand
    #print(pca_dirs.dir_mat)
    experi_strat1 = SlidingPCAStrat(unit_model, lifespan=3, pca_dirs=pca_dirs)
    experi_input1 = dict(model=unit_model,
                        strat=experi_strat1,
                        label="",
                        num_steps=NUM_STEPS)

    experi_strat2 = SlidingPCAStrat(unit_model, lifespan=3, pca_dirs=pca_dirs)
    experi_input2 = dict(model=unit_model,
                        strat=experi_strat2,
                        label="",
                        num_steps=NUM_STEPS)

    experi = PhasePlotExperiment(experi_input1, experi_input2)
    experi.execute()
    experi.plot_results(0,1, separate=False)

    Timer.generate_stats()
