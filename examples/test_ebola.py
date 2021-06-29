from kaa.reach import ReachSet
from models.ebola import Ebola

def test_Ebola():
    model = Ebola()
    mod_reach = ReachSet(model)
    mod_flow = mod_reach.computeReachSet(100)

    mod_flow.plot2DProj(0,1)
