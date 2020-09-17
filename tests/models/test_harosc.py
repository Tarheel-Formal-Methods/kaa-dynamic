from kaa.reach import ReachSet
from kaa.plotutil import Plot
from models.harosc import HarOsc

from kaa.timer import Timer

def test_HarOsc():

    model = HarOsc()
    #trajs = generate_traj(model, 10, 200)
    mod_reach = ReachSet(model)
    mod_flow = mod_reach.computeReachSet(5)

    sir_plot = Plot()
    #trajs = generate_traj(model, 10, 200)

    'Generaste the trajectories and add them to the plot.'
    sir_plot.add(mod_flow)
    sir_plot.plot(0,1)
    
    Timer.generate_stats()
