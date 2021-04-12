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

def test_sapo_skewed_compare_SIR():
    num_steps = 150
    use_supp = True
    use_pregen = False

    num_trajs = 5000
    num_steps = 70

    pca_window_size = 10
    lin_window_size = 10
    model_1 = SIR(delta=0.5)

    experi_input_1 = dict(model=model_1, #Encompass strat initilizations?
                        strat=None,
                        label="SapoSIR",
                        num_steps=num_steps-1,
                        supp_mode = use_supp,
                        pregen_mode = use_pregen,
                        num_trajs=num_trajs,
                        max_steps=num_steps)



    model = SIR_UnitBox(delta=0.5)

    pca_strat = SlidingPCAStrat(model, lifespan=pca_window_size)
    lin_strat = SlidingLinStrat(model, lifespan=lin_window_size)

    experi_input_2 = dict(model=model,
                        strat=MultiStrategy(pca_strat, lin_strat),
                        label=f"SlidingPCA Step {pca_window_size} and SlidingLin Step {lin_window_size}",
                        supp_mode = use_supp,
                        pregen_mode = use_pregen,
                        num_trajs=num_trajs,
                        num_steps=num_steps-1,
                        max_steps=num_steps)

    harosc = VolumeExperiment(experi_input_1, experi_input_2, label=f"SAPO1010COMP")
    harosc.execute(1)

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
    test_skewed_sliding_strat_comb(unit_model, 100, 5000, use_supp=True, use_pregen=False)

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
