from termcolor import colored

from kaa.reach import ReachSet
from kaa.plotutil import Plot
From models.sir import SIR_UnitBox, SIR
from models.rossler import Rossler, Rossler_UnitBox
from models.lotkavolterra import LotkaVolterra, LotkaVolterra_UnitBox
from models.quadcopter import Quadcopter, Quadcopter_UnitBox

from kaa.timer import Timer
from kaa.experiutil import generate_traj, sup_error_bounds
from kaa.temp.pca_strat import PCAStrat
from kaa.bundle import BundleMode

NUM_STEPS = 300
ITER_SPREAD = 5

def test_sir_pca_strat():

    #Compute Sapo's version.
    sir_pca = SIR_UnitBox()
    sir = SIR()
    sir_reach = ReachSet(sir)

    #sir_flow = sir_reach.computeReachSet(NUM_STEPS)
    sir_plot = Plot()
    #sir_plot.add(sir_flow)

    for i in range(3, 4):
        print(colored("Generating PCA with Iterative Step Size: {}".format(i), "white", attrs=['reverse', 'blink']))
        sir_pca_reach = ReachSet(sir_pca)
        sir_flow_pca = sir_pca_reach.computeReachSet(NUM_STEPS, tempstrat=PCAStrat(sir_pca, iter_steps=i), transmode=BundleMode.AFO)
        sir_plot.add(sir_flow_pca, "SIR_PCA_{}".format(i))
    
    sir_plot.plot(0,1,2)
    Timer.generate_stats()
    err_s = sup_error_bounds(sir_flow, sir_flow_pca, 0)
    err_i = sup_error_bounds(sir_flow, sir_flow_pca, 1)
    err_r = sup_error_bounds(sir_flow, sir_flow_pca, 2)
    print(colored("Maximum Error between PCA and Sapo along S: {}".format(err_s),'red',attrs=['blink']))
    print(colored("Maximum Error between PCA and Sapo along I: {}".format(err_i),'red',attrs=['blink']))
    print(colored("Maximum Error between PCA and Sapo along R: {}".format(err_r),'red',attrs=['blink']))

def test_rossler_pca_strat():

    #Compute Sapo's version.
    rossler_pca = Rossler_UnitBox()
    rossler = Rossler()
    rossler_reach = ReachSet(rossler)

    #rossler_flow = rossler_reach.computeReachSet(NUM_STEPS)
    rossler_plot = Plot()
    #rossler_plot.add(rossler_flow)

    for i in range(3, 4):
        print(colored("Generating PCA with Iterative Step Size: {}".format(i), "white", attrs=['reverse', 'blink']))
        rossler_pca_reach = ReachSet(rossler_pca)
        rossler_flow_pca = rossler_pca_reach.computeReachSet(NUM_STEPS, PCAStrat(rossler_pca, iter_steps=i))
        rossler_plot.add(rossler_flow_pca, "Rossler_PCA_{}".format(i))

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
