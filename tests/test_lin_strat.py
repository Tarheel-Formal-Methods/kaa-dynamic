from termcolor import colored

from kaa.reach import ReachSet
from kaa.plotutil import Plot
from models.sir import SIR_UnitBox, SIR
from models.rossler import Rossler, Rossler_UnitBox
from models.lotkavolterra import LotkaVolterra, LotkaVolterra_UnitBox
from models.quadcopter import Quadcopter, Quadcopter_UnitBox

from kaa.timer import Timer
from kaa.experiutil import generate_traj, sup_error_bounds
from kaa.temp.lin_app_strat import LinStrat
from kaa.settings import PlotSettings

NUM_STEPS = 3
ITER_SPREAD = 2

PlotSettings.save_fig = False

def test_sir_lin_strat():

    #Compute Sapo's version.
    sir_lin = SIR_UnitBox(delta=0.5)
    sir = SIR()
    #sir_reach = ReachSet(sir)

    #sir_flow = sir_reach.computeReachSet(NUM_STEPS)
    sir_plot = Plot()
    #sir_plot.add(sir_flow)

    for i in range(10,11):
        print(colored("Generating Lin_Approx with Iterative Step Size: {}".format(i), "white", attrs=['reverse', 'blink']))
        sir_lin_reach = ReachSet(sir_lin)
        sir_flow_lin = sir_lin_reach.computeReachSet(NUM_STEPS, LinStrat(sir_lin, iter_steps=i))
        sir_plot.add(sir_flow_lin, "SIR_LIN_{}".format(i))

    #sir_plot.plot(0,1,2)
    sir_plot.plot2DPhase(0,1,separate=True)
    Timer.generate_stats()

    """
    err_s = sup_error_bounds(sir_flow, sir_flow_lin, 0)
    err_i = sup_error_bounds(sir_flow, sir_flow_lin, 1)
    err_r = sup_error_bounds(sir_flow, sir_flow_lin, 2)
    print(colored("Maximum Error between Lin_Approx and Sapo along S: {}".format(err_s),'red',attrs=['blink']))
    print(colored("Maximum Error between Lin_Approx and Sapo along I: {}".format(err_i),'red',attrs=['blink']))
    print(colored("Maximum Error between Lin_Approx and Sapo along R: {}".format(err_r),'red',attrs=['blink']))
    """

def test_rossler_lin_strat():

    #Compute Sapo's version.
    rossler_lin = Rossler_UnitBox()
    rossler = Rossler()
    rossler_reach = ReachSet(rossler)

    rossler_flow = rossler_reach.computeReachSet(NUM_STEPS)
    rossler_plot = Plot()
    rossler_plot.add(rossler_flow)

    for i in range(10,11):
        print(colored("Generating LIN with Iterative Step Size: {}".format(i), "white", attrs=['reverse', 'blink']))
        rossler_lin_reach = ReachSet(rossler_lin)
        rossler_flow_lin = rossler_lin_reach.computeReachSet(NUM_STEPS, LinStrat(rossler_lin, iter_steps=i))
        rossler_plot.add(rossler_flow_lin, "Rossler_LIN_{}".format(i))

    rossler_plot.plot(0,1,2)
    Timer.generate_stats()
