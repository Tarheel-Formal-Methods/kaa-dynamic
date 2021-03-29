from models.jetengine import JetEngine_UnitBox
from kaa.temp.pca_strat import *
from kaa.temp.lin_app_strat import *
from kaa.temp.random_static_strat import *
from kaa.templates import MultiStrategy
from kaa.experiment import *
from kaa.experi_init import *

from kaa.settings import PlotSettings, KaaSettings
from kaa.timer import Timer

def test_box_JetEngine():
    num_steps = 100

    model = JetEngine_UnitBox()

    experi_input = dict(model=model,
                        strat=None,
                        label=f"JetEngine Box Reachable Set",
                        num_steps=num_steps)

    experi = PhasePlotExperiment(experi_input)
    experi.execute(0, 1, plot_border_traj=True)
    Timer.generate_stats()

def test_sliding_equal_JetEngine():
    num_steps = 250
    model = JetEngine_UnitBox(delta=0.05)
    test_equal_sliding_strat(model, num_steps)

def test_max_sliding_lin_strat_JetEngine():
    num_steps = 100
    model = JetEngine_UnitBox(delta=0.1)
    test_max_sliding_lin_strat(model, num_steps)

def test_pca_dominant_JetEngine():
    num_steps = 100
    model =  JetEngine_UnitBox()

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

def test_ran_strat_JetEngine():
    model = JetEngine_UnitBox(delta=0.1)
    test_ran_strat(model, 100, 5000, use_supp=True, use_pregen=False)

def test_skewed_sliding_strat_comb_JetEngine():
    unit_model = JetEngine_UnitBox(delta=0.1)
    test_skewed_sliding_strat_comb(unit_model, 100, 4000, use_supp=True, use_pregen=False)

def test_ran_strat_JetEngine():
    model = JetEngine_UnitBox(delta=0.1)
    test_ran_strat(model, 100, 5000, use_supp=True, use_pregen=False)

def test_sliding_pca_JetEngine():
    unit_model = JetEngine_UnitBox(delta=0.1)
    test_sliding_pca(unit_model, 20, 100, 4000, use_supp=True, use_pregen=False)

def test_sliding_lin_JetEngine():
    unit_model =JetEngine_UnitBox(delta=0.1)
    test_sliding_lin(unit_model, 20, 100, 4000, use_supp=True, use_pregen=False)
