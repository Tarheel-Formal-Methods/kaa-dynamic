from kaa.reach import ReachSet
from kaa.plotutil import Plot
from models.vanderpol import VanDerPol

from kaa.timer import Timer
def test_VDP():

    model = VanDerPol()
    mod_reach = ReachSet(model)
    mod_flow = mod_reach.computeReachSet(350)

    vdp_plot = Plot()
    vdp_plot.add(mod_flow)
    vdp_plot.plot2DPhase(0,1)

    Timer.generate_stats()
