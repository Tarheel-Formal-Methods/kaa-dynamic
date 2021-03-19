from kaa.experi_init import *

from models.LL import LL
from kaa.timer import Timer

def test_LL():
    num_steps = 200
    model = LL()

    experi_input = dict(model=model, #Encompass strat initilizations?
                        strat=None,
                        label="SapoLL",
                        num_steps=num_steps)

    harosc = ProjectionPlotExperiment(experi_input)
    harosc.execute(0, 1, 2, plot_border_traj=False)
