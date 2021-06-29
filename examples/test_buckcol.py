from models.buckcol import BuckCol_UnitBox
from kaa.experi_init import *

from kaa.timer import Timer

def test_box_BuckCol():
    num_steps = 4

    model = BuckCol_UnitBox()

    experi_input = dict(model=model,
                        strat=None,
                        label=f"BuckCol Box Reachable Set",
                        num_steps=num_steps)

    experi = PhasePlotExperiment(experi_input)
    experi.execute(0, 1, plot_border_traj=True)
    Timer.generate_stats()

def test_sliding_equal_BuckCol():
    num_steps = 60
    model =  BuckCol_UnitBox()
    test_equal_sliding_strat(model, num_steps)

def test_max_sliding_lin_strat_BuckCol():
    num_steps = 20
    model =  BuckCol_UnitBox()
    test_max_sliding_lin_strat(model, num_steps)

def test_pca_dominant_BuckCol():
    num_steps = 20
    model = BuckCol_UnitBox()

    use_supp = True
    use_pregen = False

    num_trajs = 5000

    pca_strat = SlidingPCAStrat(model, lifespan=5)
    lin_strat = SlidingLinStrat(model, lifespan=15)

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
