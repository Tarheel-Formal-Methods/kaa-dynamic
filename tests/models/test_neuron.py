from models.neuron import Neuron_UnitBox, Neuron
from kaa.temp.pca_strat import *
from kaa.temp.lin_app_strat import *
from kaa.temp.random_static_strat import *
from kaa.templates import MultiStrategy
from kaa.experiment import *
from kaa.experi_init import *

from kaa.settings import PlotSettings, KaaSettings
from kaa.timer import Timer

def test_box_Neuron():
    num_steps = 4

    model = Neuron_UnitBox()
    experi_input = dict(model=model,
                        strat=None,
                        label=f"Neuron Box Reachable Set",
                        num_steps=num_steps)

    experi = PhasePlotExperiment(experi_input)
    experi.execute(0, 1, plot_border_traj=True)
    Timer.generate_stats()

def test_init_reach_vol_vs_ran_Neuron():
    num_steps = 500
    use_supp = True
    use_pregen = False

    num_trajs = 5000

    pca_window_size = 18
    lin_window_size = 2

    inputs = []
    for inc in range(0,10,1):
        inc /= 500

        box = ((0.9-inc,1.1), (2.4-inc,2.6))
        unit_model = Neuron_UnitBox(init_box=box)
        model = Neuron(init_box=box)

        pca_strat = SlidingPCAStrat(unit_model, lifespan=pca_window_size)
        lin_strat = SlidingLinStrat(unit_model, lifespan=lin_window_size)

        experi_input_one = dict(model=unit_model,
                                strat=MultiStrategy(pca_strat, lin_strat),
                                label=f"Neuron SlidingPCA Step {pca_window_size} and SlidingLin Step {lin_window_size}",
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

def test_init_reach_vol_Neuron():
    num_steps = 200
    use_supp = True
    use_pregen = False

    num_trajs = 5000

    pca_window_size = 4
    lin_window_size = 1

    inputs_one = []
    inputs_two = []
    for inc in range(0,5,1):
        inc /= 500

        box = ((0.9-inc,1.1), (2.4-inc,2.6))
        unit_model = Neuron_UnitBox(init_box=box)
        model = Neuron(init_box=box)

        pca_strat = SlidingPCAStrat(unit_model, lifespan=pca_window_size)
        lin_strat = SlidingLinStrat(unit_model, lifespan=lin_window_size)

        experi_input_one = dict(model=unit_model,
                                strat=MultiStrategy(pca_strat, lin_strat),
                                label=f"Neuron SlidingPCA Step {pca_window_size} and SlidingLin Step {lin_window_size}",
                                supp_mode = use_supp,
                                pregen_mode = use_pregen,
                                num_trajs=num_trajs,
                                num_steps=num_steps)

        experi_input_two = dict(model=model,
                                strat=None,
                                label=f"SapoNeuron",
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

def test_pca_dominant_Neuron():
    num_steps = 500
    model =  Neuron_UnitBox()

    use_supp = True
    use_pregen = False

    num_trajs = 5000

    pca_strat = SlidingPCAStrat(model, lifespan=15)
    lin_strat = SlidingLinStrat(model, lifespan=5)

    experi_input = dict(model=model,
                        strat=MultiStrategy(pca_strat, lin_strat),
                        label=f"SlidingPCA Size 15, SlidingLin Size 5",
                        supp_mode = use_supp,
                        pregen_mode = use_pregen,
                        num_trajs=num_trajs,
                        num_steps=num_steps-1,
                        max_steps=num_steps)

    if use_supp:
        file_identifier = "(SUPP)"
    elif use_pregen:
        file_identifier = f"(PREGEN: {num_trajs})"
    else:
        file_identifier = "(RAND)"

    experi = PhasePlotExperiment(experi_input)
    experi.execute(0, 1, plot_border_traj=True)
    Timer.generate_stats()

def test_ran_strat_Neuron():
    model = Neuron_UnitBox()
    test_ran_strat(model, 500, 5000, use_supp=True, use_pregen=False)

def test_skewed_sliding_strat_comb_Neuron():
    unit_model =  Neuron_UnitBox()
    model = Neuron()
    test_skewed_sliding_strat_comb(model, 200, 5000, num_temps=5, incre=1, use_supp=True, use_pregen=False, use_sapo=model)

def test_sliding_pca_Neuron():
    model =  Neuron_UnitBox()
    test_sliding_pca(model, 20, 500, 5000, use_supp=True, use_pregen=False)

def test_sliding_lin_Neuron():
    model =  Neuron_UnitBox()
    test_sliding_lin(model, 20, 500, 5000, use_supp=True, use_pregen=False)
