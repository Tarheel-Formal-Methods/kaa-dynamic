from kaa.reach import ReachSet
from kaa.flowpipe import FlowPipePlotter
from models.lotkavolterra import LotkaVolterra

from kaa.timer import Timer
def test_LV():

    model = LotkaVolterra()
    mod_reach = ReachSet(model)
    mod_flow = mod_reach.computeReachSet(50)

    FlowPipePlotter(mod_flow).plot2DProj(0)
    FlowPipePlotter(mod_flow).plot2DProj(1)
    FlowPipePlotter(mod_flow).plot2DProj(2)

    Timer.generate_stats()
