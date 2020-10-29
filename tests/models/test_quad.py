from kaa.reach import ReachSet
from kaa.plotutil import Plot
from models.quadcopter import Quadcopter, Quadcopter_UnitBox
from kaa.temp.pca_strat import PCAStrat

from kaa.timer import Timer

def test_Quad():

    model = Quadcopter()
    mod_reach = ReachSet(model)
    mod_flow = mod_reach.computeReachSet(3)

    quad_plot = Plot()
    quad_plot.add(mod_flow)
    quad_plot.plot(2,5,13)

    Timer.generate_stats()


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
