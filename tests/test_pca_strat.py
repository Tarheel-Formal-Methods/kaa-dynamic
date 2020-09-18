from kaa.reach import ReachSet
from kaa.plotutil import Plot
from models.sir import SIR_UnitBox, SIR
from models.rossler import Rossler_UnitBox
from models.lotkavolterra import LotkaVolterra_UnitBox
from models.quadcopter import Quadcopter_UnitBox

from kaa.timer import Timer
from kaa.experiutil import generate_traj
from kaa.temp.pca_strat import PCAStrat

NUM_STEPS = 100
ITER_SPREAD = 5

def test_pca_strat():

    #Compute Sapo's version.
    sir_pca = SIR_UnitBox()
    sir = SIR()
    sir_reach = ReachSet(sir)

    sir_flow = sir_reach.computeReachSet(NUM_STEPS)
    sir_plot = Plot()
    sir_plot.add(sir_flow)

    for i in range(1,ITER_SPREAD):
        print("Generating PCA with Iterative Step Size: {}".format(2*i))
        sir_pca_reach = ReachSet(sir_pca)
        sir_flow_pca = sir_pca_reach.computeReachSet(NUM_STEPS, PCAStrat(sir_pca, iter_steps=2*i))
        sir_plot.add(sir_flow_pca, "SIR_PCA_{}".format(2*i))
    
    sir_plot.plot(0,1,2, overlap=False)
    Timer.generate_stats()
