from kaa.reach import ReachSet
from kaa.plotutil import Plot
from models.sir import SIR

from kaa.timer import Timer
from kaa.experiutil import generate_traj
from kaa.temp.simple_strat import SimpleStrat

def test_simple_strat():

    model = SIR()
    #trajs = generate_traj(model, 10, 200)
    mod_reach = ReachSet(model)
    mod_flow = mod_reach.computeReachSet(200, SimpleStrat)

    sir_plot = Plot()
    #trajs = generate_traj(model, 10, 200)

    'Generate the trajectories and add them to the plot.'
    sir_plot.add(mod_flow)

    sir_plot.plot(0,1,2)
    
    Timer.generate_stats()
