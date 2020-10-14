from kaa.reach import ReachSet
from kaa.plotutil import Plot
from models.harosc import HarOsc
from kaa.temp.pca_strat import PCAStrat
from kaa.temp.lin_app_strat import LinStrat
from kaa.settings import PlotSettings
from kaa.trajectory import Traj

from kaa.timer import Timer
from itertools import product

PlotSettings.save_fig = False

def test_HarOsc():

    model = HarOsc()
    #trajs = generate_traj(model, 10, 200)
    mod_reach = ReachSet(model)
    #mod_flow = mod_reach.computeReachSet()

    sir_plot = Plot()

    #mod_flow = mod_reach.computeReachSet(NUM_STEPS)

    SIR_PCA_ITER_STEPS = 1 #Number of steps between each recomputation of PCA Templates.
    'PCA Strategy Parameters'
    SIR_PCA_TRAJ_STEPS = 1 #Number of steps our sample trajectories should run.
    SIR_PCA_NUM_TRAJ = 100 #Number of sample trajectories we should use for the PCA routine.

    pca_strat = PCAStrat(model, traj_steps=SIR_PCA_TRAJ_STEPS, num_trajs=SIR_PCA_NUM_TRAJ, iter_steps=SIR_PCA_ITER_STEPS)
    mod_pca_flow = mod_reach.computeReachSet(3, tempstrat=pca_strat)
    #trajs = generate_traj(model, 10, 200)

    'Generaste the trajectories and add them to the plot.'
    sir_plot.add(mod_pca_flow)
    sir_plot.plot2DPhase(0,1, separate=True, plotvertices=True)
    
    Timer.generate_stats()


def test_lin_HarOsc():

    NUM_STEPS = 3

    model = HarOsc()
    #trajs = generate_traj(model, 10, 200)
    mod_reach = ReachSet(model)
    #mod_flow = mod_reach.computeReachSet()

    sir_plot = Plot()

    #mod_flow = mod_reach.computeReachSet(NUM_STEPS)

    SIR_PCA_ITER_STEPS = 1 #Number of steps between each recomputation of PCA Templates.
    'PCA Strategy Parameters'
    SIR_PCA_TRAJ_STEPS = 1 #Number of steps our sample trajectories should run.
    SIR_PCA_NUM_TRAJ = 100 #Number of sample trajectories we should use for the PCA routine.

    pca_strat = LinStrat(model, iter_steps=SIR_PCA_ITER_STEPS)
    mod_pca_flow = mod_reach.computeReachSet(NUM_STEPS, tempstrat=pca_strat)
    trajs = [Traj(model, point, steps=NUM_STEPS) for point in product([-5,-4],[0,1])]

    'Generaste the trajectories and add them to the plot.'
    sir_plot.add(mod_pca_flow, "HarOsc LINAPP")
    for t in trajs:
        sir_plot.add(t)

    sir_plot.plot2DPhase(0,1, separate=True, plotvertices=True)

    Timer.generate_stats()
