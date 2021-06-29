from models.jetengine import Lorentz_UnitBox
from kaa.experi_init import *

from kaa.timer import Timer

def test_box_Lorentz():
    num_steps = 100

    model = Lorentz_UnitBox()
    experi_input = dict(model=model,
                        strat=None,
                        label=f"Lorenz Box Reachable Set",
                        num_steps=num_steps)

    experi = PhasePlotExperiment(experi_input)
    experi.execute(0, 1, plot_border_traj=True)
    Timer.generate_stats()

def test_sliding_equal_Lorentz():
    model = Lorentz_UnitBox()
    test_equal_sliding_strat(model)
