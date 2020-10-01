from kaa.reach import ReachSet
from kaa.plotutil import Plot
from models.phos import Phosphorelay, Phosphorelay_UnitBox

from kaa.timer import Timer

def test_Phos():

    model = Phosphorelay()
    unit_model = Phosphorelay_UnitBox()
    mod_reach = ReachSet(model)
    mod_unit_reach = ReachSet(unit_model)
    unit_flow = mod_unit_reach.computeReachSet(200)
    mod_flow = mod_reach.computeReachSet(200)

    phos_plot = Plot()
    phos_plot.add(mod_flow)
    phos_plot.plot(0,1,2)

    Timer.generate_stats()
