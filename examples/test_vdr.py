from models.vanderpol import VanDerPol, VanDerPol_UnitBox
from kaa.experi_init import *
from kaa.timer import Timer
from kaa.modes import BundleTransMode


def test_sapo_VDP():
    num_steps = 70

    model = VanDerPol(delta=0.08,init_box=((0, 0.05),(1.95, 2)))
    experi_input = dict(model=model,
                        strat=None,
                        label=f"Sapo's Reachable Set",
                        num_steps=num_steps,
                        trans_mode=BundleTransMode.AFO)

    experi = PhasePlotExperiment(experi_input)
    experi.execute(0, 1)
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
                        supp_mode=use_supp,
                        pregen_mode=use_pregen,
                        num_trajs=num_trajs,
                        num_steps=num_steps-1,
                        trans_mode=BundleTransMode.AFO)

    harosc = VolumePlotExperiment(experi_input)
    harosc.execute()

def test_sliding_phase_plot_VDP():
    use_supp = True
    use_pregen = False

    num_trajs = 5000
    num_steps = 70

    model = VanDerPol_UnitBox(delta=0.08, init_box=((0, 0.05),(1.95, 2)))

    pca_window_size = 5
    lin_window_size = 0

    pca_strat = SlidingPCAStrat(model, lifespan=pca_window_size)
    lin_strat = SlidingLinStrat(model, lifespan=lin_window_size)

    experi_input = dict(model=model,
                        strat=MultiStrategy(pca_strat, lin_strat),
                        label=f"VDP SlidingPCA Step {pca_window_size} and SlidingLin Step {lin_window_size}",
                        supp_mode=use_supp,
                        pregen_mode=use_pregen,
                        num_trajs=num_trajs,
                        num_steps=num_steps)

    if use_supp:
        file_identifier = "(SUPP)"
    elif use_pregen:
        file_identifier = f"(PREGEN: {num_trajs})"
    else:
        file_identifier = "(RAND)"

    experi = PhasePlotExperiment(experi_input, separate=True)
    experi.execute(0,1)
    Timer.generate_stats()

def test_OFO_vs_AFO_phase_plot_VDP():
    use_supp = True
    use_pregen = False

    num_trajs = 5000
    num_steps = 70

    model = VanDerPol_UnitBox(delta=0.08, init_box=((0, 0.01),(1.99, 2)))

    pca_window_size = 0
    lin_window_size = 20

    pca_strat = SlidingPCAStrat(model, lifespan=pca_window_size)
    lin_strat = SlidingLinStrat(model, lifespan=lin_window_size)

    experi_input_afo = dict(model=model,
                        strat=MultiStrategy(pca_strat, lin_strat),
                        label=f"VDP AFO SlidingPCA Step {pca_window_size} and SlidingLin Step {lin_window_size}",
                        supp_mode = use_supp,
                        pregen_mode = use_pregen,
                        num_trajs=num_trajs,
                        num_steps=num_steps,
                        trans_mode=BundleTransMode.AFO)

    experi_input_ofo = dict(model=model,
                        strat=MultiStrategy(pca_strat, lin_strat),
                        label=f"VDP OFO SlidingPCA Step {pca_window_size} and SlidingLin Step {lin_window_size}",
                        supp_mode = use_supp,
                        pregen_mode = use_pregen,
                        num_trajs=num_trajs,
                        num_steps=num_steps,
                        trans_mode=BundleTransMode.OFO)

    if use_supp:
        file_identifier = "(SUPP)"
    elif use_pregen:
        file_identifier = f"(PREGEN: {num_trajs})"
    else:
        file_identifier = "(RAND)"

    experi = PhasePlotExperiment(experi_input_afo, experi_input_ofo, separate=False)
    experi.execute(0,1)
    Timer.generate_stats()

def test_init_reach_vol_vs_ran_VDP():
    num_steps = 70
    use_supp = True
    use_pregen = False

    num_trajs = 5000

    pca_window_size = 8
    lin_window_size = 12

    inputs = []
    for inc in range(5):
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

        inputs.append(experi_input_one)

    if use_supp:
        file_identifier = "(SUPP)"
    elif use_pregen:
        file_identifier = f"(PREGEN: {num_trajs})"
    else:
        file_identifier = "(RAND)"

    experi = InitReachVSRandomPlotExperiment(*inputs, num_ran_temps=pca_window_size+lin_window_size, num_trials=10)
    experi.execute()

def test_init_reach_vol_VDP():
    num_steps = 70
    use_supp = True
    use_pregen = False

    num_trajs = 5000

    pca_window_size = 3
    lin_window_size = 2

    inputs_one = []
    inputs_two = []
    for inc in range(5):
        inc /= 100

        box = ((0, 0.01+inc),(1.99 - inc, 2))

        unit_model = VanDerPol_UnitBox(delta=0.08, init_box=box)
        model = VanDerPol(delta=0.08, init_box=box)

        pca_strat = SlidingPCAStrat(unit_model, lifespan=pca_window_size)
        lin_strat = SlidingLinStrat(unit_model, lifespan=lin_window_size)

        experi_input_one = dict(model=unit_model,
                                strat=MultiStrategy(pca_strat, lin_strat),
                                label=f"VDP PCA WinSize {pca_window_size} and Lin WinSize {lin_window_size}",
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

def test_ran_diag_static_VDP():
    unit_model = VanDerPol_UnitBox(delta=0.08)
    test_ran_diag_static(unit_model, 70, 5)

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
    test_skewed_sliding_strat_comb(unit_model, 70, 4000, num_temps=5, incre=1, use_supp=True, use_pregen=False, use_sapo=model)

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
