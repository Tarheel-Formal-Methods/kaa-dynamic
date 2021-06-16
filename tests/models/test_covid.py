from kaa.reach import ReachSet
from kaa.plotutil import Plot
from models.covid import Covid, Covid_UnitBox
from kaa.temp.pca_strat import *
from kaa.temp.lin_app_strat import *
from kaa.templates import MultiStrategy

from kaa.timer import Timer
from kaa.trajectory import Traj
from kaa.experi_init import *

from kaa.bundle import BundleTransMode

def test_proj_Covid():
    num_steps = 400
    use_supp = True
    use_pregen = False

    pca_window_size = 8
    lin_window_size = 12

    num_trajs = 5000


    unit_model = Covid_UnitBox()

    pca_strat = SlidingPCAStrat(unit_model, lifespan=pca_window_size)
    lin_strat = SlidingLinStrat(unit_model, lifespan=lin_window_size)

    experi_input = dict(model=unit_model,
                        strat=MultiStrategy(pca_strat, lin_strat),
                        label=f"Covid PCA WinSize {pca_window_size} and Lin WinSize {lin_window_size}",
                        supp_mode = use_supp,
                        pregen_mode = use_pregen,
                        num_trajs=num_trajs,
                        num_steps=num_steps,
                        trans_mode = BundleTransMode.AFO)

    experi = ProjectionPlotExperiment(experi_input)
    experi.execute(0,1,2)
    Timer.generate_stats()

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
    for inc in range(4):
        inc /= 50

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

    avg_ran_vals = [3.69E-09, 1.05E-05, 0.000485901782596, 0.00638598473362]
    experi = InitReachVSRandomPlotExperiment(*inputs, num_ran_temps=pca_window_size+lin_window_size, num_trials=3, log_scale=True, precalc_vals=avg_ran_vals)
    experi.execute()

def test_init_reach_vol_Covid():
    num_steps = 200
    use_supp = True
    use_pregen = False

    num_trajs = 5000

    lin_window_size = 2
    pca_window_size = 1

    inputs_one = []
    inputs_two = []
    for inc in range(5):
        inc /= 50

        box = ((0.69 - inc, 0.70), (0.09-inc,0.1), (0.14-inc, 0.15), (0.04-inc, 0.05), (0, inc), (0, inc), (0, inc))

        unit_model = Covid_UnitBox(init_box=box)
        model = Covid(init_box=box)

        pca_strat = SlidingPCAStrat(unit_model, lifespan=pca_window_size)
        lin_strat = SlidingLinStrat(unit_model, lifespan=lin_window_size)

        experi_input_one = dict(model=unit_model,
                                strat=MultiStrategy(pca_strat, lin_strat),
                                label=f"Covid PCA WinSize {pca_window_size} and Lin WinSize {lin_window_size}",
                                supp_mode = use_supp,
                                pregen_mode = use_pregen,
                                num_trajs=num_trajs,
                                num_steps=num_steps)

        experi_input_two = dict(model=model,
                                strat=None,
                                label=f"SapoCovid",
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

    experi = InitReachPlotExperiment(*inputs, log_scale=True)
    experi.execute()

def test_covid_data_plot():
    experi = CovidDataPlotExperiment(earliest_date='06/21/20',
                                     latest_date='08/22/20',
                                     beta_interval=(0.10, 0.12),
                                     gamma_interval=(0.078, 0.082),
                                     eta=0.0015)
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
    test_skewed_sliding_strat_comb(unit_model, 200, 5000, num_temps=3, incre=1, use_supp=True, use_pregen=False, use_sapo=model)

def test_sliding_strat_comb_Covid():
    model = Covid_UnitBox()
    test_sliding_strat_comb(model, 150, 5000, use_supp=True, use_pregen=False)
