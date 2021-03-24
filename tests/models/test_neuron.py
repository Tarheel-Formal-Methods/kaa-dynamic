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
    num_steps = 300

    model = Neuron_UnitBox(delta=0.08)

    experi_input = dict(model=model,
                        strat=None,
                        label=f"Neuron Box Reachable Set",
                        num_steps=num_steps)

    experi = PhasePlotExperiment(experi_input)
    experi.execute(0, 1, plot_border_traj=True)
    Timer.generate_stats()
