from models.chem import Chem_UnitBox
from kaa.experiment import *
from kaa.experi_init import *
from kaa.timer import Timer
from kaa.bundle import BundleTransMode

def test_arch_Chem():
    use_supp = True
    use_pregen = False

    num_trajs = 5000

    num_steps = 40
    delta = 0.005

    model_one = Chem_UnitBox(100, 1000, delta=delta)
    #model_two = CoupledVDP_UnitBox(delta=delta, init_box=init_box , mu=2)

    pca_window_size = 5
    lin_window_size = 5

    pca_strat_one = SlidingPCAStrat(model_one, lifespan=pca_window_size)
    lin_strat_one = SlidingLinStrat(model_one, lifespan=lin_window_size)

    experi_input_one = dict(model=model_one,
                        strat=MultiStrategy(pca_strat_one, lin_strat_one),
                        label=f"ROBE21",
                        num_steps=num_steps,
                        supp_mode = use_supp,
                        pregen_mode = use_pregen,
                        num_trajs=num_trajs,
                        trans_mode=BundleTransMode.AFO)

    experi = ProjectionPlotExperiment(experi_input_one, plot_total_width=True)
    experi.execute()
