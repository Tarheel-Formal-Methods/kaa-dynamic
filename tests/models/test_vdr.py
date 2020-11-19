from kaa.reach import ReachSet
from kaa.plotutil import Plot
from kaa.trajectory import Traj
from models.vanderpol import VanDerPol, VanDerPol_UnitBox

from kaa.temp.pca_strat import PCAStrat, DelayedPCAStrat
from kaa.temp.lin_app_strat import LinStrat
from kaa.templates import MultiStrategy
from kaa.experiment import PhasePlotExperiment

from kaa.settings import PlotSettings
from kaa.timer import Timer

PlotSettings.save_fig = False


def test_VDP():
    NUM_STEPS = 70

    model = VanDerPol(delta=0.08)
    vdp_sapo = PhasePlotExperiment(model)
    vdp_sapo.execute(NUM_STEPS)

    vdp_sapo.plot_results(0,1)
    Timer.generate_stats()


def test_pca_VDP():

    NUM_STEPS = 4

    unit_model = VanDerPol_UnitBox()

    VDP_PCA_ITER_STEPS = 1 #Number of steps between each recomputation of PCA Templates.
    'PCA Strategy Parameters'
    VDP_PCA_TRAJ_STEPS = 5 #Number of steps our sample trajectories should run.
    VDP_PCA_NUM_TRAJ = 100 #Number of sample trajectories we should use for the PCA routine.

    pca_strat = PCAStrat(unit_model, traj_steps=VDP_PCA_TRAJ_STEPS, num_trajs=VDP_PCA_NUM_TRAJ, iter_steps=VDP_PCA_ITER_STEPS)
    vdp_pca = PhasePlotExperiment(unit_model, pca_strat)
    vdp_pca.execute(NUM_STEPS)
    vdp_pca.plot_results(0,1)

    Timer.generate_stats()


def test_lin_VDP():

    NUM_STEPS = 50
    VDP_LIN_ITER_STEPS = 1 #Number of steps between each recomputation of LinApp Templates.

    unit_model = VanDerPol_UnitBox(delta=0.02)
    unit_mod_reach = ReachSet(unit_model)

    lin_strat = LinStrat(unit_model, iter_steps=VDP_LIN_ITER_STEPS)
    mod_lin_flow = unit_mod_reach.computeReachSet(NUM_STEPS, tempstrat=lin_strat)

    #points = [[0,1.97], [0.01, 1.97], [0.01,2], [0,2], [0.005,1.97], [0.005,2], [0,1.97],  [0,1.985], [0.01,1.985]]
    #trajs = [Traj(unit_model, point, NUM_STEPS) for point in points]

    vdp_plot = Plot()
    vdp_plot.add(mod_lin_flow, "VDP LinAPP")

    'Add trajectories'
    #for traj in trajs:
    #    vdp_plot.add(traj)

    vdp_plot.plot2DPhase(0,1, separate=True)

    Timer.generate_stats()

def test_pca_lin_VDP():

    NUM_STEPS = 70
    VDP_LIN_ITER_STEPS = 1 #Number of steps between each recomputation of LinApp Templates.
    VDP_PCA_ITER_STEPS = 1 #Number of steps between each recomputation of PCA Templates.
    'PCA Strategy Parameters'
    VDP_PCA_TRAJ_STEPS = 5 #Number of steps our sample trajectories should run.
    VDP_PCA_NUM_TRAJ = 200 #Number of sample trajectories we should use for the PCA routine.
    VDP_PCA_DELAY = 5

    unit_model = VanDerPol_UnitBox(delta=0.08)

    lin_strat = MultiStrategy(LinStrat(unit_model, iter_steps=VDP_LIN_ITER_STEPS), \
                              PCAStrat(unit_model, traj_steps=VDP_PCA_TRAJ_STEPS, num_trajs=VDP_PCA_NUM_TRAJ, iter_steps=VDP_PCA_ITER_STEPS), \
                              PCAStrat(unit_model, traj_steps=VDP_PCA_TRAJ_STEPS, num_trajs=VDP_PCA_NUM_TRAJ, iter_steps=VDP_PCA_ITER_STEPS+VDP_PCA_DELAY),
                              PCAStrat(unit_model, traj_steps=VDP_PCA_TRAJ_STEPS, num_trajs=VDP_PCA_NUM_TRAJ, iter_steps=VDP_PCA_ITER_STEPS+2*VDP_PCA_DELAY))

    vdp_pca = PhasePlotExperiment(unit_model, lin_strat)
    vdp_pca.execute(NUM_STEPS)
    vdp_pca.plot_results(0,1)

    Timer.generate_stats()



def test_delayed_pca_VDP():

    NUM_STEPS = 70
    VDP_LIN_ITER_STEPS = 1 #Number of steps between each recomputation of LinApp Templates.
    VDP_PCA_ITER_STEPS = 1 #Number of steps between each recomputation of PCA Templates.
    'PCA Strategy Parameters'
    VDP_PCA_TRAJ_STEPS = 5 #Number of steps our sample trajectories should run.
    VDP_PCA_NUM_TRAJ = 200 #Number of sample trajectories we should use for the PCA routine.
    #
    VDP_PCA_LIFE_SPAN = 3

    unit_model = VanDerPol_UnitBox(delta=0.08)

    lin_strat = MultiStrategy(LinStrat(unit_model, iter_steps=VDP_LIN_ITER_STEPS), \
                              DelayedPCAStrat(unit_model, traj_steps=VDP_PCA_TRAJ_STEPS, num_trajs=VDP_PCA_NUM_TRAJ, life_span=VDP_PCA_LIFE_SPAN))

    vdp_pca = PhasePlotExperiment(unit_model, lin_strat)
    vdp_pca.execute(NUM_STEPS)
    vdp_pca.plot_results(0,1)

    Timer.generate_stats()
