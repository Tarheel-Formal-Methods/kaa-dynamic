from termcolor import colored

from kaa.reach import ReachSet
from kaa.plotutil import Plot
from models.sir import SIR_UnitBox, SIR
from models.rossler import Rossler, Rossler_UnitBox
from models.lotkavolterra import LotkaVolterra, LotkaVolterra_UnitBox
from models.quadcopter import Quadcopter, Quadcopter_UnitBox

from kaa.timer import Timer
from kaa.temp.pca_strat import PCAStrat
from kaa.temp.lin_app_strat import LinStrat
from kaa.bundle import BundleMode

NUM_STEPS = 300

def pca_lin_comp():

    sir = SIR_UnitBox()
    sir_reach = ReachSet(sir)

    sir_flow = sir_reach.computeReachSet(NUM_STEPS)
    sir_plot = Plot()
    sir_plot.add(sir_flow)

    sir_lin_flow = sir_reach.computeReachSet(NUM_STEPS, LinStrat(sir, iter_steps=5))
    sir_pca_flow = sir_reach.computeReachSet(NUM_STEPS, PCAStrat(sir, iter_steps=5))

    sir_plot.add(sir_lin_flow, "SIR_LIN")
    sir_plot.add(sir_pca_flow, "SIR_PCA")

    sir_plot.plot(0,1,2)
    Timer.generate_stats()
