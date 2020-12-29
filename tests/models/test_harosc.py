from models.harosc import HarOsc
from kaa.temp.pca_strat import *
from kaa.temp.lin_app_strat import *
from kaa.templates import MultiStrategy
from kaa.experiment import Experiment, PhasePlotExperiment

from kaa.settings import PlotSettings
from itertools import product

PlotSettings.save_fig = False

def test_pca_HarOsc():

    NUM_STEPS = 4

    model = HarOsc()
    #trajs = generate_traj(model, 10, 200)
    mod_reach = ReachSet(model)
    #mod_flow = mod_reach.computeReachSet()

    sir_plot = Plot()

    SIR_PCA_ITER_STEPS = 1 #Number of steps between each recomputation of PCA Templates.
    'PCA Strategy Parameters'
    SIR_PCA_TRAJ_STEPS = 2 #Number of steps our sample trajectories should run.
    SIR_PCA_NUM_TRAJ = 100 #Number of sample trajectories we should use for the PCA routine.

    pca_strat = PCAStrat(model, traj_steps=SIR_PCA_TRAJ_STEPS, num_trajs=SIR_PCA_NUM_TRAJ, iter_steps=SIR_PCA_ITER_STEPS)
    mod_pca_flow = mod_reach.computeReachSet(NUM_STEPS, tempstrat=[pca_strat], transmode=BundleMode.AFO)
    #trajs = generate_traj(model, 10, 200)

    'Generaste the trajectories and add them to the plot.'
    sir_plot.add(mod_pca_flow, "HarOsc PCA")
    sir_plot.plot2DPhase(0,1, separate=True, plotvertices=True)
    
    Timer.generate_stats()


def test_lin_HarOsc():

    NUM_STEPS = 5

    model = HarOsc()
    #trajs = generate_traj(model, 10, 200)
    mod_reach = ReachSet(model)
    #mod_flow = mod_reach.computeReachSet()

    sir_plot = Plot()

    #mod_flow = mod_reach.computeReachSet(NUM_STEPS)

    SIR_LIN_ITER_STEPS = 1 #Number of steps between each recomputation of PCA Templates.

    lin_strat = LinStrat(model, iter_steps=SIR_LIN_ITER_STEPS)
    mod_lin_flow = mod_reach.computeReachSet(NUM_STEPS, tempstrat=lin_strat, transmode=BundleMode.AFO)
    trajs = [Traj(model, point, steps=NUM_STEPS) for point in product([-5,-4],[0,1])]

    'Generaste the trajectories and add them to the plot.'
    sir_plot.add(mod_lin_flow, "HarOsc LINAPP")
    for t in trajs:
        sir_plot.add(t)

    sir_plot.plot2DPhase(0,1, separate=True, plotvertices=True)

    Timer.generate_stats()


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
    #print(pca_dirs.dir_mat)
    #print(len(pca_dirs.dir_mat))
    experi.plot_results(0,1, separate=False)

    Timer.generate_stats()
