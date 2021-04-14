from itertools import product

from models.vanderpol import VanDerPol, VanDerPol_UnitBox
from kaa.temp.pca_strat import *
from kaa.temp.lin_app_strat import *
from kaa.temp.random_static_strat import *
from kaa.templates import MultiStrategy
from kaa.experiment import *
from kaa.experi_init import *

from kaa.settings import PlotSettings, KaaSettings
from kaa.timer import Timer

PlotSettings.save_fig = False
def test_sapo_VDP():
    num_steps = 30

    model = VanDerPol(delta=0.08)

    experi_input = dict(model=model,
                        strat=None,
                        label=f"Sapo's Reachable Set",
                        num_steps=num_steps)

    experi = PhasePlotExperiment(experi_input)
    experi.execute(0, 1, plot_border_traj=True)
    Timer.generate_stats()

def test_sapo_vol_VDP():
    use_supp = True
    use_pregen = False

    num_trajs = 5000
    num_steps = 70

    model = VanDerPol(delta=0.08)
    experi_input = dict(model=model, #Encompass strat initilizations?
                        strat=None,
                        label="SapoVDP",
                        supp_mode = use_supp,
                        pregen_mode = use_pregen,
                        num_trajs=num_trajs,
                        num_steps=num_steps-1,
                        max_steps=num_steps)

    harosc = VolumePlotExperiment(experi_input)
    harosc.execute()

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

def test_init_reach_vol_VDP():
    num_steps = 30
    use_supp = True
    use_pregen = False

    num_trajs = 5000

    pca_window_size = 8
    lin_window_size = 12

    inputs_one = []
    inputs_two = []
    for inc in range(3):
        inc /= 100

        box = ((0, 0.01+inc),(1.99 - inc, 2))

        unit_model = VanDerPol_UnitBox(delta=0.08, init_box=box)
        model = VanDerPol(delta=0.08, init_box=box)

        pca_strat = SlidingPCAStrat(unit_model, lifespan=pca_window_size)
        lin_strat = SlidingLinStrat(unit_model, lifespan=lin_window_size)

        experi_input_one = dict(model=unit_model,
                                strat=MultiStrategy(pca_strat, lin_strat),
                                label=f"VDP SlidingPCA Step {pca_window_size} and SlidingLin Step {lin_window_size}",
                                supp_mode = use_supp,
                                pregen_mode = use_pregen,
                                num_trajs=num_trajs,
                                num_steps=num_steps)

        experi_input_two = dict(model=model,
                                strat=None,
                                label=f"SapoVDP",
                                supp_mode = use_supp,
                                pregen_mode = use_pregen,
                                num_trajs=num_trajs,
                                num_steps=num_steps)


        inputs_one.append(experi_input_one)
        inputs_two.append(experi_input_two)


    inputs = inputs_one + inputs_two

    if use_supp:
        file_identifier = "(SUPP)"
    elif use_pregen:
        file_identifier = f"(PREGEN: {num_trajs})"
    else:
        file_identifier = "(RAND)"

    experi = InitReachPlotExperiment(*inputs)
    experi.execute()

def test_ani_lin_comp_VDP():
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

def test_vol_comp_VDP():
    unit_model = VanDerPol_UnitBox(delta=0.08)
    model = VanDerPol(delta=0.08)
    test_vol_comp(unit_model, 70, 4000, use_supp=True, use_pregen=False, use_sapo=None)

def test_sliding_equal_VDP():
    model = VanDerPol_UnitBox(delta=0.08)
    test_equal_sliding_strat(model)

def test_ran_strat_VDP():
    model = VanDerPol_UnitBox(delta=0.08)
    test_ran_strat(model, 70, 5000, use_supp=True, use_pregen=False)

def test_strat_comb_VDP():
    unit_model = VanDerPol_UnitBox(delta=0.08)
    test_strat_comb(unit_model, (1,3,5), 70, 4000, use_supp=True, use_pregen=False)

def test_skewed_sliding_strat_comb_VDP():
    unit_model = VanDerPol_UnitBox(delta=0.08)
    model = VanDerPol(delta=0.08)
    test_skewed_sliding_strat_comb(unit_model, 70, 4000, use_supp=True, use_pregen=False, use_sapo=model)

def test_sliding_strat_comb_VDP():
    unit_model = VanDerPol_UnitBox(delta=0.08)
    test_sliding_strat_comb(unit_model, 70, 4000, use_supp=True, use_pregen=False)

def test_one_one_strat_pca_VDP():
    unit_model = VanDerPol_UnitBox(delta=0.08)
    test_one_one_strat_pca(unit_model, 70, 4000)

def test_one_one_strat_lin_VDP():
    unit_model = VanDerPol_UnitBox(delta=0.08)
    test_one_one_strat_lin(unit_model, 70, 4000)

def test_sliding_pca_VDP():
    unit_model = VanDerPol_UnitBox(delta=0.08)
    test_sliding_pca(unit_model, 20, 1, 4000, use_supp=True, use_pregen=False)

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
