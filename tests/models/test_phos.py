from kaa.reach import ReachSet
from kaa.plotutil import Plot
from models.phos import Phosphorelay

import kaa.benchmark as Timer

def test_Phos():

    model = Phosphorelay()
    mod_reach = ReachSet(model)
    mod_flow = mod_reach.computeReachSet(200)

    phos_plot = Plot()
    phos_plot.add(mod_flow)
    phos_plot.plot(0,1,2)

    Timer.generate_stats()
