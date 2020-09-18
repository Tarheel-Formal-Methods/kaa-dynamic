from termcolor import colored

from kaa.reach import ReachSet
from kaa.plotutil import Plot
from models.sir import SIR_UnitBox, SIR
from models.rossler import Rossler, Rossler_UnitBox
from models.lotkavolterra import LotkaVolterra, LotkaVolterra_UnitBox
from models.quadcopter import Quadcopter, Quadcopter_UnitBox

from kaa.timer import Timer
from kaa.experiutil import generate_traj
from kaa.temp.pca_strat import PCAStrat

NUM_STEPS = 100
ITER_SPREAD = 2

def test_sir_pca_strat():

    #Compute Sapo's version.
    sir_pca = SIR_UnitBox()
    sir = SIR()
    sir_reach = ReachSet(sir)

    sir_flow = sir_reach.computeReachSet(NUM_STEPS)
    sir_plot = Plot()
    sir_plot.add(sir_flow)

    for i in range(1,ITER_SPREAD):
        print(colored("Generating PCA with Iterative Step Size: {}".format(2*i), "white", attrs=['reverse', 'blink']))
        sir_pca_reach = ReachSet(sir_pca)
        sir_flow_pca = sir_pca_reach.computeReachSet(NUM_STEPS, PCAStrat(sir_pca, iter_steps=2*i))
        sir_plot.add(sir_flow_pca, "SIR_PCA_{}".format(2*i))
    
    sir_plot.plot(0,1,2)
    Timer.generate_stats()

def test_rossler_pca_strat():

    #Compute Sapo's version.
    rossler_pca = Rossler_UnitBox()
    rossler = Rossler()
    rossler_reach = ReachSet(rossler)

    rossler_flow = rossler_reach.computeReachSet(NUM_STEPS)
    rossler_plot = Plot()
    rossler_plot.add(rossler_flow)

    for i in range(1,ITER_SPREAD):
        print(colored("Generating PCA with Iterative Step Size: {}".format(2*i), "white", attrs=['reverse', 'blink']))
        rossler_pca_reach = ReachSet(rossler_pca)
        rossler_flow_pca = rossler_pca_reach.computeReachSet(NUM_STEPS, PCAStrat(rossler_pca, iter_steps=2*i))
        rossler_plot.add(rossler_flow_pca, "Rossler_PCA_{}".format(2*i))

    rossler_plot.plot(0,1,2)
    Timer.generate_stats()

def test_lv_pca_strat():

    #Compute Sapo's version.
    lv_pca = LotkaVolterra_UnitBox()
    lv = LotkaVolterra()
    lv_reach = ReachSet(lv)

    lv_flow = lv_reach.computeReachSet(NUM_STEPS)
    lv_plot = Plot()
    lv_plot.add(lv_flow)

    for i in range(1,ITER_SPREAD):
        print(colored("Generating PCA with Iterative Step Size: {}".format(2*i), "white", attrs=['reverse','blink']))
        lv_pca_reach = ReachSet(lv_pca)
        lv_flow_pca = lv_pca_reach.computeReachSet(NUM_STEPS, PCAStrat(lv_pca, iter_steps=2*i))
        lv_plot.add(lv_flow_pca, "LV_PCA_{}".format(2*i))

    lv_plot.plot(0,1,2)
    Timer.generate_stats()

def test_quad_pca_strat():

    #Compute Sapo's version.
    quad_pca = Quadcopter_UnitBox()
    quad = Quadcopter()
    quad_reach = ReachSet(quad)

    quad_flow = quad_reach.computeReachSet(NUM_STEPS)
    quad_plot = Plot()
    quad_plot.add(quad_flow)

    for i in range(1,ITER_SPREAD):
        print(colored("Generating PCA with Iterative Step Size: {}".format(2*i), "white", attrs=['reverse','blink']))
        quad_pca_reach = ReachSet(quad_pca)
        quad_flow_pca = quad_pca_reach.computeReachSet(NUM_STEPS, PCAStrat(quad_pca, iter_steps=2*i))
        quad_plot.add(quad_flow_pca, "QUAD_PCA_{}".format(2*i))

    quad_plot.plot(2,5,13)
    Timer.generate_stats()
