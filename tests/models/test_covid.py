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

def test_init_reach_vol_vs_ran_Covid():
    num_steps = 200
    use_supp = True
    use_pregen = False

    num_trajs = 5000

    pca_window_size = 8
    lin_window_size = 12

    inputs = []
    for inc in range(5):
        inc /= 100

        box = ((0.69 -inc, 0.70), (0.09-inc,0.1), (0.14-inc, 0.15), (0.04-inc, 0.05), (0.00099, 0.001), (0.00099, 0.001), (0.00099, 0.001))

        unit_model = Covid_UnitBox(init_box=box)

        pca_strat = SlidingPCAStrat(unit_model, lifespan=pca_window_size)
        lin_strat = SlidingLinStrat(unit_model, lifespan=lin_window_size)

        experi_input_one = dict(model=unit_model,
                                strat=MultiStrategy(pca_strat, lin_strat),
                                label=f"COVID SlidingPCA Step {pca_window_size} and SlidingLin Step {lin_window_size}",
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


def test_equal_sliding_strat_Covid():
    model = Covid_UnitBox()
    test_equal_sliding_strat(model, 100, vars=(3,4,6))

def test_ran_strat_Covid():
    model = Covid_UnitBox()
    test_ran_strat(model, 200, 5000, use_supp = True, use_pregen = False)

def test_skewed_sliding_strat_comb_Covid():
    unit_model = Covid_UnitBox()
    model = Covid()
    test_skewed_sliding_strat_comb(unit_model, 200, 5000, num_temps=3, incre=1, use_supp=True, use_pregen=False, use_sapo=)

def test_sliding_strat_comb_Covid():
    model = Covid_UnitBox()
    test_sliding_strat_comb(model, 150, 5000, use_supp=True, use_pregen=False)
