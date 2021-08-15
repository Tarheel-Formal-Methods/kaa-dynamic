from models.invertpend import InvertPend, InvertPend_UnitBox
from kaa.temp.pca_strat import *
from kaa.temp.lin_app_strat import *
from kaa.templates import MultiStrategy
from kaa.experi_init import *
from kaa.experi_lib import *
from kaa.modes import BundleTransMode

from settings import PlotSettings
PlotSettings.save_fig = False

def test_InvertPend():
    num_steps = 80

    init_box = ((0.25,0.3),(0.25, 0.3))
    model = InvertPend(init_box=init_box)

    pca_window_size = 10
    lin_window_size = 10

    pca_strat_one = SlidingPCAStrat(model, lifespan=pca_window_size)
    lin_strat_one = SlidingLinStrat(model, lifespan=lin_window_size)

    experi_input = dict(model=model,
                        strat=MultiStrategy(pca_strat_one, lin_strat_one),
                        label="InvertPend",
                        supp_mode=True,
                        pregen_mode=True,
                        num_steps=num_steps,
                        trans_mode=BundleTransMode.AFO)

    harosc = PhasePlotExperiment(experi_input)
    harosc.execute(0,1)

def test_skewed_sliding_strat_comb_InvertPend():
    unit_model = InvertPend_UnitBox()
    model = InvertPend()
    test_skewed_sliding_strat_comb(unit_model, 9, 4000,
                                  num_temps=5,
                                  incre=1,
                                  use_supp=True,
                                  use_pregen=False,
                                  use_sapo=model,
                                  mode="Data")

def test_ran_strat_InvertPend():
    model = InvertPend_UnitBox()
    test_ran_strat(model, 80, 5000, use_supp=True, use_pregen=False)

def test_vol_comp_InvertPend():
    unit_model = InvertPend_UnitBox()
    model = InvertPend()
    test_vol_comp(unit_model, 80, 4000, use_supp=True, use_pregen=False, use_sapo=model)
