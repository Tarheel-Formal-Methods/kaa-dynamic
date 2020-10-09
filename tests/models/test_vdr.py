from kaa.reach import ReachSet
from kaa.plotutil import Plot
from models.vanderpol import VanDerPol, VanDerPol_UnitBox

from kaa.temp.pca_strat import PCAStrat
from kaa.temp.lin_app_strat import LinStrat

from kaa.settings import PlotSettings
from kaa.timer import Timer

PlotSettings.save_fig = False

def test_VDP():

    NUM_STEPS = 100

    model = VanDerPol()
    unit_model = VanDerPol_UnitBox()

    mod_reach = ReachSet(model)
    unit_mod_reach = ReachSet(unit_model)
    mod_flow = mod_reach.computeReachSet(NUM_STEPS)

    VDP_PCA_ITER_STEPS = 100 #Number of steps between each recomputation of PCA Templates.
    'PCA Strategy Parameters'
    VDP_PCA_TRAJ_STEPS = 5 #Number of steps our sample trajectories should run.
    VDP_PCA_NUM_TRAJ = 100 #Number of sample trajectories we should use for the PCA routine.

    pca_strat = PCAStrat(unit_model, traj_steps=VDP_PCA_TRAJ_STEPS, num_trajs=VDP_PCA_NUM_TRAJ, iter_steps=VDP_PCA_ITER_STEPS)
    mod_pca_flow = unit_mod_reach.computeReachSet(NUM_STEPS, tempstrat=pca_strat)

    vdp_plot = Plot()
    vdp_plot.add(mod_flow, "VDP SAPO")
    vdp_plot.add(mod_pca_flow, "VDP PCA")
    vdp_plot.plot2DPhase(0,1, separate=False, plotvertices=True)

    Timer.generate_stats()


def test_lin_VDP():

    NUM_STEPS = 100

    model = VanDerPol()
    unit_model = VanDerPol_UnitBox()

    mod_reach = ReachSet(model)
    unit_mod_reach = ReachSet(unit_model)
    mod_flow = mod_reach.computeReachSet(NUM_STEPS)

    VDP_PCA_ITER_STEPS = 10 #Number of steps between each recomputation of PCA Templates.
    'PCA Strategy Parameters'
    VDP_PCA_TRAJ_STEPS = 5 #Number of steps our sample trajectories should run.
    VDP_PCA_NUM_TRAJ = 100 #Number of sample trajectories we should use for the PCA routine.

    pca_strat = LinStrat(unit_model, iter_steps=VDP_PCA_ITER_STEPS)
    mod_pca_flow = unit_mod_reach.computeReachSet(NUM_STEPS, tempstrat=pca_strat)

    vdp_plot = Plot()
    vdp_plot.add(mod_flow)
    vdp_plot.add(mod_pca_flow)
    vdp_plot.plot2DPhase(0,2, separate=False, plotvertices=True)

    Timer.generate_stats()
