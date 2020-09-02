from kaa.reach import ReachSet
from kaa.plotutil import Plot
from models.lotkavolterra import LotkaVolterra

from kaa.timer import Timer
def test_LV():

    model = LotkaVolterra()
    mod_reach = ReachSet(model)
    mod_flow = mod_reach.computeReachSet(200)

    plot = Plot()
    plot.add(mod_flow)
    plot.plot(0,1,2)

    Timer.generate_stats()
