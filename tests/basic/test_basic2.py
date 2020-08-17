from kaa.reach import ReachSet
from kaa.timer import Timer
from kaa.plotutil import Plot
from models.basic.basic2 import Basic2


def test_basic2():

    basic_mod = Basic2()
    basic_reach = ReachSet(basic_mod)
    flowpipe = basic_reach.computeReachSet(300)
    
    basic_plot = Plot()
    basic_plot.add(flowpipe)
    basic_plot.plot(0)

    Timer.generate_stats()
