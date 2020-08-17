from kaa.reach import ReachSet
from kaa.plotutil import Plot
from models.basic.basic import Basic

def test_basic():

    basic_mod = Basic()
    basic_reach = ReachSet(basic_mod)
    flowpipe = basic_reach.computeReachSet(300)

    basic_plot = Plot()
    basic_plot.add_flowpipe(flowpipe)
    basic_plot.plot(0)
