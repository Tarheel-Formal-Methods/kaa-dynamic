from kaa.reach import ReachSet
from kaa.plotutil import Plot
from models.phos import Phosphorelay, Phosphorelay_UnitBox
from kaa.templates import MultiStrategy
from kaa.temp.lin_app_strat import LinStrat
from kaa.temp.pca_strat import PCAStrat, DelayedPCAStrat
from kaa.templates import MultiStrategy

from kaa.timer import Timer
from kaa.trajectory import Traj

def test_Phos():

    model = Phosphorelay()
    #unit_model = Phosphorelay_UnitBox()
    mod_reach = ReachSet(model)
    #mod_unit_reach = ReachSet(unit_model)
    #unit_flow = mod_unit_reach.computeReachSet(200)
    mod_flow = mod_reach.computeReachSet(30)

    phos_plot = Plot()
    phos_plot.add(mod_flow)
    phos_plot.plot2DPhase(0,1, separate=False, plotvertices=True)

    Timer.generate_stats()


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
