from models.invertpend import InvertPend
from kaa.temp.pca_strat import *
from kaa.temp.lin_app_strat import *
from kaa.templates import MultiStrategy
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