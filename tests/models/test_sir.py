from kaa.reach import ReachSet
from kaa.plotutil import Plot
from models.sir import SIR, SIR_UnitBox
from kaa.temp.pca_strat import *
from kaa.temp.lin_app_strat import *
from kaa.templates import MultiStrategy

from kaa.timer import Timer
from kaa.trajectory import Traj
from kaa.experiutil import *

from kaa.bundle import BundleMode

def test_SIR():

    model = SIR(delta=0.5)
    model_unit = SIR_UnitBox()
    #trajs = generate_traj(model, 10, 200)
    mod_reach = ReachSet(model)
    mod_unit_reach = ReachSet(model_unit)
    mod_flow = mod_reach.computeReachSet(70)
    #mod_unit_flow = mod_unit_reach.computeReachSet(300)

    sir_plot = Plot()
    #trajs = generate_traj(model, 10, 200)

    'Generaste the trajectories and add them to the plot.'
    #for traj in trajs:
    #    sir_plot.add(traj)
    sir_plot.add(mod_flow)
    #sir_plot.add(mod_unit_flow)
    sir_plot.plot2DPhase(0,1)
    
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
    sir_plot.plot2DPhase(1,2,separate=False, plotvertices=True)
    sir_plot.plot2DPhase(0,2,separate=False, plotvertices=True)

    Timer.generate_stats()


def test_strat_comb_sir():
    unit_model = SIR_UnitBox(delta=0.5)
    test_strat_comb(unit_model, (1,3,5), 70)
