from kaa.reach import ReachSet
from kaa.plotutils import Plot
from models.quadcopter import Quadcopter

from kaa.timer import Timer

def test_Quad():

    model = Quadcopter()
    mod_reach = ReachSet(model)
    mod_flow = mod_reach.computeReachSet(10)

    quad_plot = Plot()
    quad_plot(mod_flow)
    quad_plot.plot(2,5,13)

    Timer.generate_stats()
