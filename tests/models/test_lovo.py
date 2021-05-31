from models.lovo import LOVO_UnitBox
from kaa.experiment import *
from kaa.experi_init import *
from kaa.timer import Timer
from kaa.bundle import BundleTransMode

def test_arch_LOVO21():
    use_supp = True
    use_pregen = False

    num_trajs = 5000

    num_steps = 120
    delta = 0.03

    model_one = LOVO_UnitBox(delta=delta)

    pca_window_size = 0
    lin_window_size = 20

    pca_strat_one = SlidingPCAStrat(model_one, lifespan=pca_window_size)
    lin_strat_one = SlidingLinStrat(model_one, lifespan=lin_window_size)

    experi_input_one = dict(model=model_one,
                            strat=MultiStrategy(pca_strat_one, lin_strat_one),
                            label=f"LOVO21",
                            num_steps=num_steps,
                            supp_mode = use_supp,
                            pregen_mode = use_pregen,
                            num_trajs=num_trajs,
                            trans_mode=BundleTransMode.AFO)

    experi = PhasePlotExperiment(experi_input_one)
    experi.execute(0, 1)
