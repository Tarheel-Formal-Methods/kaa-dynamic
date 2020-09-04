from kaa.reach import ReachSet
from kaa.plotutil import Plot
from models.sir import SIR_UnitBox, SIR
from models.rossler import Rossler_UnitBox
from models.lotkavolterra import LotkaVolterra_UnitBox
from models.quadcopter import Quadcopter_UnitBox

from kaa.timer import Timer
from kaa.experiutil import generate_traj
from kaa.temp.pca_strat import PCAStrat

def test_pca_strat():

    sir_pca = SIR_UnitBox()
    sir = SIR()
    sir_pca_reach = ReachSet(sir_pca)
    sir_reach = ReachSet(sir)

    sir_flow_pca = sir_pca_reach.computeReachSet(200, PCAStrat)
    print("Generating SIR")
    sir_flow = sir_reach.computeReachSet(200)
    sir_plot = Plot()

    sir_plot.add(sir_flow)
    sir_plot.add(sir_flow_pca, ``"SIR_PCA")
    sir_plot.plot(0,1,2)

    Timer.generate_stats()
