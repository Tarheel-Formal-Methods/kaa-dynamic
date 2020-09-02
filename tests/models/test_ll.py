from kaa.reach import ReachSet
from kaa.plotutil import Plot
from models.LL import LL

from kaa.timer import Timer

def test_LL():

    model = LL()
    mod_reach = ReachSet(model)
    mod_flow = mod_reach.computeReachSet(150)

    ll_plot = Plot()
    ll_plot.add(mod_flow)
    ll_plot.plot(0)

    Timer.generate_stats()
