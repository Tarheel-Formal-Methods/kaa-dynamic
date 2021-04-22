from kaa.reach import ReachSet
from kaa.plotutil import Plot
from models.sir import SIR, SIR_UnitBox
from kaa.temp.pca_strat import *
from kaa.temp.lin_app_strat import *
from kaa.templates import MultiStrategy

from kaa.timer import Timer
from kaa.trajectory import Traj
from kaa.experi_init import *

from kaa.bundle import BundleMode

def test_sapo_SIR():
    num_steps = 150
    model = SIR(delta=0.5)

    experi_input = dict(model=model, #Encompass strat initilizations?
                        strat=None,
                        label="SapoSIR",
                        supp_mode = False,
                        pregen_mode = False,
                        num_trajs=5000,
                        num_steps=num_steps-1,
                        max_steps=num_steps)

    harosc = ProjectionPlotExperiment(experi_input)
    harosc.execute(0,1,2)
    Timer.generate_stats()

def test_sapo_vol_SIR():
    use_supp = True
    use_pregen = False

    num_trajs = 5000
    num_steps = 150

    model = SIR(delta=0.5)
    experi_input = dict(model=model, #Encompass strat initilizations?
                        strat=None,
                        label="SapoSIR",
                        supp_mode = use_supp,
                        pregen_mode = use_pregen,
                        num_trajs=num_trajs,
                        num_steps=num_steps-1,
                        max_steps=num_steps)

    harosc = VolumeExperiment(experi_input)
    harosc.execute(1)

def test_init_reach_vol_vs_ran_SIR():
    num_steps = 150
    use_supp = True
    use_pregen = False

    num_trajs = 5000

    pca_window_size = 8
    lin_window_size = 12

    inputs = []
    for inc in range(0,5,1):
        inc /= 500

        box = ((0.79-inc,0.8), (0.19-inc,0.2), (0.00099, 0.001))

        unit_model = SIR_UnitBox(init_box=box)
        model = SIR(init_box=box)

        pca_strat = SlidingPCAStrat(unit_model, lifespan=pca_window_size)
        lin_strat = SlidingLinStrat(unit_model, lifespan=lin_window_size)

        experi_input_one = dict(model=unit_model,
                                strat=MultiStrategy(pca_strat, lin_strat),
                                label=f"SIR SlidingPCA Step {pca_window_size} and SlidingLin Step {lin_window_size}",
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

def test_init_reach_vol_vs_sapo_SIR():
    num_steps = 150
    use_supp = True
    use_pregen = False

    num_trajs = 5000

    pca_window_size = 8
    lin_window_size = 12

    inputs_one = []
    inputs_two = []
    for inc in range(0,5,1):
        inc /= 500

        box = ((0.79-inc,0.8), (0.19-inc,0.2), (0.00099, 0.001))

        unit_model = SIR_UnitBox(init_box=box)
        model = SIR(init_box=box)

        pca_strat = SlidingPCAStrat(unit_model, lifespan=pca_window_size)
        lin_strat = SlidingLinStrat(unit_model, lifespan=lin_window_size)

        experi_input_one = dict(model=unit_model,
                                strat=MultiStrategy(pca_strat, lin_strat),
                                label=f"{model.name} SlidingPCA Step {pca_window_size} and SlidingLin Step {lin_window_size}",
                                supp_mode = use_supp,
                                pregen_mode = use_pregen,
                                num_trajs=num_trajs,
                                num_steps=num_steps)

        experi_input_two = dict(model=model,
                                strat=None,
                                label=f"SapoSIR",
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

def test_vol_comp_SIR():
    unit_model = SIR_UnitBox()
    model = SIR()
    test_vol_comp(unit_model, 100, 4000, use_supp=True, use_pregen=False, use_sapo=None)

def test_ran_strat_SIR():
    model = SIR_UnitBox()
    test_ran_strat(model, 150, 5000, use_supp = True, use_pregen = False)

def test_strat_comb_SIR():
    model = SIR_UnitBox()
    test_strat_comb(model, (1,3,5), 150, -1)

def test_skewed_sliding_strat_comb_SIR():
    unit_model = SIR_UnitBox()
    model = SIR()
    test_skewed_sliding_strat_comb(unit_model, 150, 5000, use_supp=True, use_pregen=False, use_sapo=model)

def test_sliding_strat_comb_SIR():
    model = SIR_UnitBox()
    test_sliding_strat_comb(model, 150, 5000, use_supp=True, use_pregen=False)

def test_one_one_strat_pca_SIR():
    model = SIR_UnitBox()
    test_one_one_strat_pca(model, 150)

def test_one_one_strat_lin_SIR():
    model = SIR_UnitBox()
    test_one_one_strat_lin(model, 150)

def test_sliding_pca_SIR():
    model = SIR_UnitBox()
    test_sliding_pca(model, 20, 150, 5000, use_supp=True, use_pregen=False)

def test_sliding_lin_SIR():
    model = SIR_UnitBox()
    test_sliding_lin(model, 20, 150, 5000, use_supp=True, use_pregen=False)

def gen_save_dirs_SIR():
    model = SIR_UnitBox(d)
    gen_save_dirs(model, 150)

def find_pca_variation_SIR():
    unit_model = SIR_UnitBox()
    find_pca_variation(unit_model, 150, max_num_trajs=6000, label="PCADevSIR")

def find_lin_variation_SIR():
    unit_model = SIR_UnitBox()
    find_lin_variation(unit_model, 150, max_num_trajs=6000, label="LinDevSIR")
