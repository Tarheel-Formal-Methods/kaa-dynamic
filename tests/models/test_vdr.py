from itertools import product

from models.vanderpol import VanDerPol, VanDerPol_UnitBox
from kaa.temp.pca_strat import *
from kaa.temp.lin_app_strat import *
from kaa.templates import MultiStrategy
from kaa.experiment import *
from kaa.experiutil import *

from kaa.settings import PlotSettings, KaaSettings
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
    pca_strat = PCAStrat(unit_model, traj_steps=VDP_PCA_TRAJ_STEPS, num_trajs=VDP_PCA_NUM_TRAJ, iter_steps=VDP_PCA_ITER_STEPS, pca_dirs=pca_dirs)

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

    NUM_STEPS = 10
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

    experi_input1 = dict(model=model,
                   strat = None,
                   label="VDP Sapo",
                   num_steps=NUM_STEPS)

    experi_input2 = dict(model=unit_model,
                   strat=lin_strat,
                   label="",
                   num_steps=NUM_STEPS)

    vdp_pca = PhasePlotExperiment(experi_input2)
    vdp_pca.execute()
    vdp_pca.plot_results(0,1)

    Timer.generate_stats()

def plot_sliding_lin_VDP():
    NUM_STEPS = 70
    NUM_TRAJS = 5000
    unit_model = VanDerPol_UnitBox(delta=0.08)

    ran_pca_strat = SlidingLinStrat(unit_model, lifespan=20)
    ran_experi_input = dict(model=unit_model, #Encompass strat initilizations?
                            strat=ran_pca_strat,
                            label="Lin Pre-gen 10 Lifespan",
                            num_steps=NUM_STEPS,
                            max_steps=70,
                            num_trajs=NUM_TRAJS,
                            supp_mode=False,
                            pregen_mode=True)

    vdp_pca = PhasePlotExperiment(ran_experi_input)
    vdp_pca.execute(0,1)

def test_ani_pca_comp_VDP():
    NUM_STEPS = 5
    NUM_TRAJS = 5000
    unit_model = VanDerPol_UnitBox(delta=0.08)

    ran_pca_strat = SlidingPCAStrat(unit_model, lifespan=20)
    ran_experi_input = dict(model=unit_model, #Encompass strat initilizations?
                            strat=ran_pca_strat,
                            label="PCA Pre-gen",
                            num_steps=NUM_STEPS,
                            max_steps=70,
                            num_trajs=NUM_TRAJS,
                            supp_mode=False,
                            pregen_mode=True)

    supp_pca_strat = SlidingPCAStrat(unit_model, lifespan=20)
    supp_experi_input = dict(model=unit_model,
                            strat=supp_pca_strat,
                            label="PCA Supp Points",
                            num_steps=NUM_STEPS,
                            max_steps=70,
                            num_trajs=NUM_TRAJS,
                            supp_mode=True,
                            pregen_mode=False)

    vdp_pca = CompAniExperiment(ran_experi_input, supp_experi_input)
    vdp_pca.execute(0, 1, 1, "VDPPCAComp", plot_pts=[False, True])

    Timer.generate_stats()

def test_ani_lin_ wd wdcomp_VDP():
    NUM_STEPS = 40
    NUM_TRAJS = 1000
    unit_model = VanDerPol_UnitBox(delta=0.08)

    ran_pca_strat = SlidingLinStrat(unit_model, lifespan=20, num_trajs=NUM_TRAJS)
    ran_experi_input = dict(model=unit_model,
                            strat=ran_pca_strat,
                            label="LinApp Pre-gen",
                            num_steps=NUM_STEPS,
                            max_steps=70,
                            num_trajs=NUM_TRAJS,
                            supp_mode=False,
                            pregen_mode=True)

    supp_pca_strat = SlidingLinStrat(unit_model, lifespan=20)
    supp_experi_input = dict(model=unit_model,
                            strat=supp_pca_strat,
                            label="LinApp Supp Points",
                            num_steps=NUM_STEPS,
                            max_steps=70,
                            num_trajs=NUM_TRAJS,
                            supp_mode=True,
                            pregen_mode=False)
                             
    vdp_pca = CompAniExperiment(ran_experi_input, supp_experi_input)
    vdp_pca.execute(0, 1, 0, "VDPLinComp", plot_pts=[False, True])

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

    lin_dirs = GeneratedLinDirs(unit_model, NUM_STEPS+1)
    lin_1 = LinStrat(unit_model, iter_steps=VDP_LIN_ITER_STEPS)
    lin_2 = LinStrat(unit_model, iter_steps=VDP_LIN_ITER_STEPS+VDP_LIN_DELAY)

    lin_strat = MultiStrategy(lin_1, lin_2)
    inputs =  dict(model=unit_model,
                   strat=lin_strat,
                   label="",
                   num_steps=NUM_STEPS)

    vdp_pca = Animation(inputs)
    vdp_pca.execute()
    vdp_pca.animate(0,1, lin_1, lin_2)
    Timer.generate_stats()

def test_strat_comb_VDP():
    unit_model = VanDerPol_UnitBox(delta=0.08)
    test_strat_comb(unit_model, (1,3,5), 70, 4000)

def test_one_one_strat_pca_VDP():
    unit_model = VanDerPol_UnitBox(delta=0.08)
    test_one_one_strat_pca(unit_model, 70, 4000)

def test_one_one_strat_lin_VDP():
    unit_model = VanDerPol_UnitBox(delta=0.08)
    test_one_one_strat_lin(unit_model, 70, 4000)

def test_sliding_pca_VDP():
    unit_model = VanDerPol_UnitBox(delta=0.08)
    test_sliding_pca(unit_model, 20, 70, 4000)

def test_sliding_lin_VDP():
    unit_model = VanDerPol_UnitBox(delta=0.08)
    test_sliding_lin(unit_model, 20, 70, 4000)

def test_strat_comb_stdev_VDP():
    unit_model = VanDerPol_UnitBox(delta=0.08)
    test_comb_stdev_reduction(unit_model, 70)

def gen_save_dirs_VDP():
    unit_model = VanDerPol_UnitBox(delta=0.08)
    gen_save_dirs(unit_model, 70, max_num_trajs=2000, num_trials=10)

def find_pca_variation_VDP():
    unit_model = VanDerPol_UnitBox(delta=0.08)
    find_pca_variation(unit_model, 70, max_num_trajs=8000)

def find_lin_variation_VDP():
    unit_model = VanDerPol_UnitBox(delta=0.08)
    find_lin_variation(unit_model, 70, max_num_trajs=8000)
