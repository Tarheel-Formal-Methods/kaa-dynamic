from kaa.reach import ReachSet
from kaa.plotutil import Plot
from models.rossler import Rossler, Rossler_UnitBox

from kaa.timer import Timer

def test_Rossler():

    model = Rossler()
    mod_reach = ReachSet(model)
    mod_flow = mod_reach.computeReachSet(300)

    rossler_plot = Plot()
    rossler_plot.add(mod_flow)
    rossler_plot.plot(0,1,2)

    Timer.generate_stats()
