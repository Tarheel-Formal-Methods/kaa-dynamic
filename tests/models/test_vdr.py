from kaa.reach import ReachSet
from models.vanderpol import VanDerPol

def test_VDP():

    model = VanDerPol()
    mod_reach = ReachSet(model)
    mod_flow = mod_reach.computeReachSet(100)

    mod_flow.plot2DPhase(0,1)
