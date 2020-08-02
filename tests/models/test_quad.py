from kaa.reach import ReachSet
from models.quadcopter import Quadcopter

from kaa.timer import Timer

def test_Quad():

    model = Quadcopter()
    mod_reach = ReachSet(model)
    mod_flow = mod_reach.computeReachSet(10)

    mod_flow.plot2DProj(2)
    mod_flow.plot2DProj(5)
    mod_flow.plot2DProj(13)

    Timer.generate_stats()
