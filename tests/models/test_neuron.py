from models.neuron import Neuron_UnitBox
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

def test_sliding_equal_Neuron():
    num_steps = 500
    model =  Neuron_UnitBox()
    test_equal_sliding_strat(model, num_steps)

def test_max_sliding_lin_strat_Neuron():
    num_steps = 500
    model =  Neuron_UnitBox()
    test_max_sliding_lin_strat(model, num_steps)

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
    model =  Neuron_UnitBox()
    test_skewed_sliding_strat_comb(model, 500, 5000, use_supp=True, use_pregen=False)

def test_sliding_pca_Neuron():
    model =  Neuron_UnitBox()
    test_sliding_pca(model, 20, 500, 5000, use_supp=True, use_pregen=False)

def test_sliding_lin_Neuron():
    model =  Neuron_UnitBox()
    test_sliding_lin(model, 20, 500, 5000, use_supp=True, use_pregen=False)
