from kaa.experi_init import *

from models.LL import LL, LL_UnitBox
from kaa.timer import Timer
from kaa.bundle import BundleTransMode

def test_LL():
    num_steps = 200
    model = LL()

    experi_input = dict(model=model, #Encompass strat initilizations?
                        strat=None,
                        label="SapoLL",
                        num_steps=num_steps)

    harosc = ProjectionPlotExperiment(experi_input)
    harosc.execute(0, 1, 2, plot_border_traj=False)

def test_arch_LL():
    use_supp = True
    use_pregen = False

    num_trajs = 5000

    num_steps = 120
    delta = 0.1

    init_box_one = ((1.19, 1.21), (1.04, 1.06), (1.49, 1.51), (2.39, 2.41), (0.99, 1.01), (0.09, 0.11), (0.44, 0.46))

    model_one = LL_UnitBox(delta=delta, init_box=init_box_one)

    pca_window_size = 5
    lin_window_size = 5

    pca_strat_one = SlidingPCAStrat(model_one, lifespan=pca_window_size)
    lin_strat_one = SlidingLinStrat(model_one, lifespan=lin_window_size)

    experi_input_one = dict(model=model_one,
                            strat=MultiStrategy(pca_strat_one, lin_strat_one),
                            label=f"LALO21",
                            num_steps=num_steps,
                            supp_mode = use_supp,
                            pregen_mode = use_pregen,
                            num_trajs=num_trajs,
                            trans_mode=BundleTransMode.AFO)

    experi = ProjectionPlotExperiment(experi_input_one)
    experi.execute(3)
