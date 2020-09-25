from termcolor import colored

from kaa.reach import ReachSet
from kaa.plotutil import Plot
from models.sir import SIR_UnitBox, SIR
from models.rossler import Rossler, Rossler_UnitBox
from models.lotkavolterra import LotkaVolterra, LotkaVolterra_UnitBox
from models.quadcopter import Quadcopter, Quadcopter_UnitBox

from kaa.timer import Timer
from kaa.experiutil import generate_traj, sup_error_bounds
from kaa.temp.cut_strat import CutStrat
from kaa.bundle import BundleMode

NUM_STEPS = 300

def test_sir_cut_strat():

    #Compute Sapo's version.
    sir_cut = SIR()
    sir = SIR_UnitBox()
    sir_reach = ReachSet(sir)

    sir_flow = sir_reach.computeReachSet(NUM_STEPS)
    sir_plot = Plot()
    sir_plot.add(sir_flow)

    sir_cut_reach = ReachSet(sir_cut)
    sir_cut_flow = sir_cut_reach.computeReachSet(NUM_STEPS, CutStrat(sir_cut))

    sir_sapo_flow = sir_cut_reach.computeReachSet(NUM_STEPS)


    sir_plot.add(sir_cut_flow, "SIR_CUT_150")
    sir_plot.add(sir_sapo_flow, "SIR_SAPO")

    sir_plot.plot(0,1,2)
    Timer.generate_stats()
