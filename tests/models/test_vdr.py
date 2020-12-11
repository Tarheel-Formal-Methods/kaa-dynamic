from itertools import product

from models.vanderpol import VanDerPol, VanDerPol_UnitBox

from kaa.temp.pca_strat import PCAStrat, DelayedPCAStrat, GeneratedPCADirs
from kaa.temp.lin_app_strat import LinStrat
from kaa.templates import MultiStrategy
from kaa.experiment import *

from kaa.settings import PlotSettings
from kaa.timer import Timer

PlotSettings.save_fig = False


def test_VDP():
    NUM_STEPS = 1

    model = VanDerPol(delta=0.08)

    vdp_sapo = PhasePlotExperiment([ExperimentInput(model, label="VDP Sapo")])
    vdp_sapo.execute(NUM_STEPS)

    vdp_sapo.plot_results(0,1)
    Timer.generate_stats()


def test_pca_VDP():

    NUM_STEPS = 3
    
    model = VanDerPol(delta=0.08)
    unit_model = VanDerPol_UnitBox(delta=0.08)

    VDP_PCA_ITER_STEPS = 1 #Number of steps between each recomputation of PCA Templates.
    'PCA Strategy Parameters'
    VDP_PCA_TRAJ_STEPS = 5 #Number of steps our sample trajectories should run.
    VDP_PCA_NUM_TRAJ = 50 #Number of sample trajectories we should use for the PCA routine.
    VDP_PCA_DELAY = 5

    pca_dirs = GeneratedPCADirs(model, VDP_PCA_NUM_TRAJ, NUM_STEPS)
    pca_strat = MultiStrategy(PCAStrat(unit_model, traj_steps=VDP_PCA_TRAJ_STEPS, num_trajs=VDP_PCA_NUM_TRAJ, iter_steps=VDP_PCA_ITER_STEPS, pca_dirs=pca_dirs), \
                              PCAStrat(unit_model, traj_steps=VDP_PCA_TRAJ_STEPS, num_trajs=VDP_PCA_NUM_TRAJ, iter_steps=VDP_PCA_ITER_STEPS+VDP_PCA_DELAY, pca_dirs=pca_dirs), \
                              PCAStrat(unit_model, traj_steps=VDP_PCA_TRAJ_STEPS, num_trajs=VDP_PCA_NUM_TRAJ, iter_steps=VDP_PCA_ITER_STEPS+2*VDP_PCA_DELAY, pca_dirs=pca_dirs))

    inputs = [ExperimentInput(model, label="VDP Sapo"), ExperimentInput(unit_model, strat=pca_strat, label="VDP Kaa PCA")]
    
    vdp_pca = PhasePlotExperiment(inputs)
    vdp_pca.execute(NUM_STEPS)
    vdp_pca.plot_results(0,1)

    Timer.generate_stats()

def test_lin_VDP():

    NUM_STEPS = 70
    VDP_LIN_ITER_STEPS = 1 #Number of steps between each recomputation of LinApp Templates.
    VDP_LIN_DELAY = 1

    model = VanDerPol(delta=0.08)
    unit_model = VanDerPol_UnitBox(delta=0.08)

    lin_strat = MultiStrategy(LinStrat(unit_model, iter_steps=VDP_LIN_ITER_STEPS), \
                              LinStrat(unit_model, iter_steps=VDP_LIN_ITER_STEPS+VDP_LIN_DELAY), \
                              LinStrat(unit_model, iter_steps=VDP_LIN_ITER_STEPS+2*VDP_LIN_DELAY))

    inputs = [ExperimentInput(model, label="VDP Sapo"), ExperimentInput(unit_model, strat=lin_strat, label="VDP Kaa Lin")]

    vdp_pca = PhasePlotExperiment(inputs)
    vdp_pca.execute(NUM_STEPS)
    vdp_pca.plot_results(0,1)

    Timer.generate_stats()

def test_pca_lin_VDP():

    NUM_STEPS = 70
    VDP_LIN_ITER_STEPS = 1 #Number of steps between each recomputation of LinApp Templates.
    VDP_PCA_ITER_STEPS = 1 #Number of steps between each recomputation of PCA Templates.
    'PCA Strategy Parameters'
    VDP_PCA_TRAJ_STEPS = 2 #Number of steps our sample trajectories should run.
    VDP_PCA_NUM_TRAJ = 200 #Number of sample trajectories we should use for the PCA routine.
    VDP_LIN_DELAY = 2
    VDP_PCA_DELAY = 5

    model = VanDerPol(delta=0.08)
    unit_model = VanDerPol_UnitBox(delta=0.08)

    lin_strat = MultiStrategy(LinStrat(unit_model, iter_steps=VDP_LIN_ITER_STEPS), \
                              LinStrat(unit_model, iter_steps=VDP_LIN_ITER_STEPS+VDP_LIN_DELAY), \
                              PCAStrat(unit_model, traj_steps=VDP_PCA_TRAJ_STEPS, num_trajs=VDP_PCA_NUM_TRAJ, iter_steps=VDP_PCA_ITER_STEPS), \
                              PCAStrat(unit_model, traj_steps=VDP_PCA_TRAJ_STEPS, num_trajs=VDP_PCA_NUM_TRAJ, iter_steps=VDP_PCA_ITER_STEPS+VDP_PCA_DELAY))

    inputs = [ExperimentInput(model, label="VDP Sapo"), ExperimentInput(unit_model, strat=lin_strat, label="VDP Kaa")]
    vdp_pca = PhasePlotExperiment(inputs)
    vdp_pca.execute(NUM_STEPS)
    vdp_pca.plot_results(0,1)

    Timer.generate_stats()

def test_delayed_pca_VDP():

    NUM_STEPS = 70
    VDP_LIN_ITER_STEPS = 1 #Number of steps between each recomputation of LinApp Templates.
    VDP_PCA_ITER_STEPS = 1 #Number of steps between each recomputation of PCA Templates.
    'PCA Strategy Parameters'
    VDP_PCA_TRAJ_STEPS = 2 #Number of steps our sample trajectories should run.
    VDP_PCA_NUM_TRAJ = 200 #Number of sample trajectories we should use for the PCA routine.
    #
    VDP_PCA_LIFE_SPAN = 3

    model = VanDerPol(delta=0.08)
    unit_model = VanDerPol_UnitBox(delta=0.08)

    lin_strat = MultiStrategy(LinStrat(unit_model, iter_steps=VDP_LIN_ITER_STEPS), \
                              DelayedPCAStrat(unit_model, traj_steps=VDP_PCA_TRAJ_STEPS, num_trajs=VDP_PCA_NUM_TRAJ, lifespan=VDP_PCA_LIFE_SPAN))

    inputs = [ExperimentInput(model, label="VDP Sapo"), ExperimentInput(unit_model, strat=lin_strat, label="VDP Kaa Delay")]
    vdp_pca = PhasePlotExperiment(inputs)
    vdp_pca.execute(NUM_STEPS)
    vdp_pca.plot_results(0,1)

    Timer.generate_stats()

def test_ani_pca_VDP():

    NUM_STEPS = 70

#model = VanDerPol(delta=0.08)
    unit_model = VanDerPol_UnitBox(delta=0.08)

    VDP_PCA_ITER_STEPS = 1 #Number of steps between each recomputation of PCA Templates.
    'PCA Strategy Parameters'
    VDP_PCA_TRAJ_STEPS = 5 #Number of steps our sample trajectories should run.
    VDP_PCA_NUM_TRAJ = 100 #Number of sample trajectories we should use for the PCA routine.
    VDP_PCA_DELAY = 5

    pca_1 = PCAStrat(unit_model, traj_steps=VDP_PCA_TRAJ_STEPS, num_trajs=VDP_PCA_NUM_TRAJ, iter_steps=VDP_PCA_ITER_STEPS)
    pca_2 = PCAStrat(unit_model, traj_steps=VDP_PCA_TRAJ_STEPS, num_trajs=VDP_PCA_NUM_TRAJ, iter_steps=VDP_PCA_ITER_STEPS+VDP_PCA_DELAY)
    pca_3 = PCAStrat(unit_model, traj_steps=VDP_PCA_TRAJ_STEPS, num_trajs=VDP_PCA_NUM_TRAJ, iter_steps=VDP_PCA_ITER_STEPS+2*VDP_PCA_DELAY)

    pca_strat = MultiStrategy(pca_1, pca_2, pca_3)

    experi_input = ExperimentInput(unit_model, strat=pca_strat, label="VDP Kaa PCA")

    vdp_pca = Animation(experi_input)
    vdp_pca.execute(NUM_STEPS)
    vdp_pca.animate(0,1, pca_1)
    vdp_pca.animate(0,1, pca_2)
    vdp_pca.animate(0,1, pca_3)

    Timer.generate_stats()


def test_ani_pca_lin_VDP():

    NUM_STEPS = 70
    VDP_LIN_ITER_STEPS = 1 #Number of steps between each recomputation of LinApp Templates.
    VDP_PCA_ITER_STEPS = 1 #Number of steps between each recomputation of PCA Templates.
    'PCA Strategy Parameters'
    VDP_PCA_TRAJ_STEPS = 2 #Number of steps our sample trajectories should run.
    VDP_PCA_NUM_TRAJ = 200 #Number of sample trajectories we should use for the PCA routine.
    VDP_LIN_DELAY = 2
    VDP_PCA_DELAY = 5

    model = VanDerPol(delta=0.08)
    unit_model = VanDerPol_UnitBox(delta=0.08)

    lin_1 = LinStrat(unit_model, iter_steps=VDP_LIN_ITER_STEPS)
    lin_2 = LinStrat(unit_model, iter_steps=VDP_LIN_ITER_STEPS+VDP_LIN_DELAY)
    pca_1 = PCAStrat(unit_model, traj_steps=VDP_PCA_TRAJ_STEPS, num_trajs=VDP_PCA_NUM_TRAJ, iter_steps=VDP_PCA_ITER_STEPS)
    pca_2 = PCAStrat(unit_model, traj_steps=VDP_PCA_TRAJ_STEPS, num_trajs=VDP_PCA_NUM_TRAJ, iter_steps=VDP_PCA_ITER_STEPS+VDP_PCA_DELAY)

    lin_strat = MultiStrategy(lin_1, lin_2, pca_1, pca_2)

    inputs = ExperimentInput(unit_model, strat=lin_strat, label="VDP Kaa")
    vdp_pca = Animation(inputs)
    vdp_pca.execute(NUM_STEPS)
    vdp_pca.animate(0,1, lin_1, lin_2)
    #vdp_pca.animate(0,1, lin_2)
    #vdp_pca.animate(0,1, pca_2)

    Timer.generate_stats()

def test_strat_comb_VDP():
    unit_model = VanDerPol_UnitBox(delta=0.08)

    NUM_STEPS = 70
    VDP_PCA_NUM_TRAJ = 300 #Number of sample trajectories we should use for the PCA routine.
    VDP_LIN_ITER_STEPS = 1 #Number of steps between each recomputation of LinApp Templates.
    VDP_PCA_ITER_STEPS = 1 #Number of steps between each recomputation of PCA Templates.
    'PCA Strategy Parameters'
    VDP_LIN_ITER_STEPS = 1
    VDP_PCA_DELAY = 2

    pca_dirs = GeneratedPCADirs(unit_model, VDP_PCA_NUM_TRAJ, NUM_STEPS+1) #Way to deduce lengeth beforehand
    pca_iter_steps = [VDP_PCA_ITER_STEPS, VDP_PCA_ITER_STEPS+VDP_PCA_DELAY]
    lin_iter_steps = [VDP_LIN_ITER_STEPS, VDP_PCA_ITER_STEPS+VDP_PCA_DELAY]

    batch1 = ExperimentBatch()
    for pca_step in range(1,6,2): #unit_model tossed around too many times.
        experi_strat = PCAStrat(unit_model, iter_steps=pca_step, pca_dirs=pca_dirs)
        experi_input = dict(model=unit_model,
                            strat=experi_strat,
                            label="",
                            num_steps=NUM_STEPS)
        
        experi = Experiment(experi_input, label=f"PCA with Box for {pca_step} steps")
        batch1.add_experi(experi)

    batch2 = ExperimentBatch()
    for pca_step in range(1,6,2): #unit_model tossed around too many times.
        experi_strat = LinStrat(unit_model, iter_steps=pca_step)
        experi_input = dict(model=unit_model,
                            strat=experi_strat,
                            label="",
                            num_steps=NUM_STEPS)
        experi = Experiment(experi_input, label=f"LinApp with Box for {pca_step} steps")
        batch2.add_experi(experi)

    exec_plot_vol_results(batch1, batch2)


def test_pca_life_VDP():
    unit_model = VanDerPol_UnitBox(delta=0.08)

    NUM_STEPS = 70
    VDP_PCA_NUM_TRAJ = 300 #Number of sample trajectories we should use for the PCA routine.
    VDP_LIN_ITER_STEPS = 1 #Number of steps between each recomputation of LinApp Templates.
    VDP_PCA_ITER_STEPS = 1 #Number of steps between each recomputation of PCA Templates.
    'PCA Strategy Parameters'
    VDP_LIN_ITER_STEPS = 2
    VDP_PCA_DELAY = 5

    LIFE_MAX = 30
    LIFE_INCREMENT = 5

    pca_dirs = GeneratedPCADirs(unit_model, VDP_PCA_NUM_TRAJ, NUM_STEPS+1) #Way to deduce lengeth beforehand
    experi_list = []

    for lifespan in range(LIFE_MAX, 0, -LIFE_INCREMENT): #unit_model tossed around too many times.
        strat = DelayedPCAStrat(unit_model, lifespan=lifespan, pca_dirs=pca_dirs)
        
        experi_input = ExperimentInput(unit_model, strat)
        experi = Experiment([experi_input])
        experi_list.append(experi)

    vol_experi_batch = VolumeExperimentBatch(experi_list)
    vol_experi_batch.execute(NUM_STEPS)
    vol_experi_batch.plot_results()
