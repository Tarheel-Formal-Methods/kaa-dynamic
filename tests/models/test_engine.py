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
    num_steps = 100
    model = JetEngine_UnitBox()
    test_equal_sliding_strats(model, num_steps)

def test_max_sliding_lin_strat_JetEngine():
    num_steps = 100
    model = JetEngine_UnitBox(delta=0.1)
    test_max_sliding_lin_strat(model, num_steps)
