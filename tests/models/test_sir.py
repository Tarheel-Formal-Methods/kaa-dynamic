from kaa.reach import ReachSet
from kaa.plotutil import Plot
from models.sir import SIR, SIR_UnitBox

from kaa.timer import Timer
from kaa.experiutil import generate_traj

from kaa.bundle import BundleMode

def test_SIR():

    model = SIR()
    model_unit = SIR_UnitBox()
    #trajs = generate_traj(model, 10, 200)
    mod_reach = ReachSet(model)
    mod_unit_reach = ReachSet(model_unit)
    mod_flow = mod_reach.computeReachSet(300)
    mod_unit_flow = mod_unit_reach.computeReachSet(300)

    sir_plot = Plot()
    #trajs = generate_traj(model, 10, 200)

    'Generaste the trajectories and add them to the plot.'
    #for traj in trajs:
    #    sir_plot.add(traj)
    sir_plot.add(mod_flow)
    sir_plot.add(mod_unit_flowq)
    sir_plot.plot(0,1,2)
    
    Timer.generate_stats()
