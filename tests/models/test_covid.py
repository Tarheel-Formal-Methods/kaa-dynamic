from kaa.reach import ReachSet
from kaa.plotutil import Plot
from models.covid import Covid, Covid_UnitBox
from kaa.temp.pca_strat import *
from kaa.temp.lin_app_strat import *
from kaa.templates import MultiStrategy

from kaa.timer import Timer
from kaa.trajectory import Traj
from kaa.experi_init import *

from kaa.bundle import BundleMode

def test_sapo_Covid():
    num_steps = 150
    model = Covid(delta=0.5)

    experi_input = dict(model=model, #Encompass strat initilizations?
                        strat=None,
                        label="SapoCovid",
                        supp_mode = False,
                        pregen_mode = False,
                        num_trajs=5000,
                        num_steps=num_steps-1,
                        max_steps=num_steps)

    harosc = ProjectionPlotExperiment(experi_input)
    harosc.execute(0,1,2)
    Timer.generate_stats()

def test_sapo_vol_Covid():
    use_supp = True
    use_pregen = False

    num_trajs = 5000
    num_steps = 150

    model = Covid(delta=0.5)
    experi_input = dict(model=model, #Encompass strat initilizations?
                        strat=None,
                        label="SapoCovid",
                        supp_mode = use_supp,
                        pregen_mode = use_pregen,
                        num_trajs=num_trajs,
                        num_steps=num_steps-1,
                        max_steps=num_steps)

    harosc = VolumeExperiment(experi_input)
    harosc.execute(1)

def test_sapo_skewed_compare_Covid():
    num_steps = 100
    use_supp = True
    use_pregen = False

    num_trajs = 5000
    num_steps = 70

    pca_window_size = 10
    lin_window_size = 10
    model_1 = Covid(delta=0.5)

    experi_input_1 = dict(model=model_1, #Encompass strat initilizations?
                        strat=None,
                        label="SapoCovid",
                        num_steps=num_steps-1,
                        supp_mode = use_supp,
                        pregen_mode = use_pregen,
                        num_trajs=num_trajs,
                        max_steps=num_steps)



    model = Covid_UnitBox(delta=0.5)

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

    harosc = Experiment(experi_input_1, experi_input_2, label=f"SAPO1010COMP")
    harosc.execute(1)

def test_equal_sliding_strat_Covid():
    model = Covid_UnitBox()
    test_equal_sliding_strat(model, 100, vars=(3,4,6))

def test_ran_strat_Covid():
    model = Covid_UnitBox()
    test_ran_strat(model, 200, 5000, use_supp = True, use_pregen = False)

def test_skewed_sliding_strat_comb_Covid():
    unit_model = Covid_UnitBox()
    model = Covid()
    test_skewed_sliding_strat_comb(unit_model, 200, 5000, use_supp=True, use_pregen=False, use_sapo=model)

def test_sliding_strat_comb_Covid():
    model = Covid_UnitBox()
    test_sliding_strat_comb(model, 150, 5000, use_supp=True, use_pregen=False)
