from kaa.reach import ReachSet
from kaa.plotutil import Plot
from settings import PlotSettings

from models.sir import SIR, SIR_UnitBox

PlotSettings.save_fig = False

def test_phase_sir():

    sir_mod = SIR_UnitBox(delta=0.5)
    sir_reach = ReachSet(sir_mod)

    flowpipe = sir_reach.computeReachSet(50)
    plot = Plot()
    plot.add(flowpipe)
    plot.plot2DPhase(1,2)
