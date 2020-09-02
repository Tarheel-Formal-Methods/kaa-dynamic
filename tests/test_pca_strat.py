from kaa.reach import ReachSet
from kaa.plotutil import Plot
from models.sir import SIR_PCA
from models.rossler import Rossler_PCA
from models.lotkavolterra import LotkaVolterra_PCA
from models.quadcopter import Quadcopter_PCA



from kaa.timer import Timer
from kaa.experiutil import generate_traj
from kaa.temp.pca_strat import PCAStrat

def test_pca_strat():

    #pca_model = SIR_PCA()
    #stan_model = SIR()
    ross_pca = Quadcopter_PCA()

    #trajs = generate_traj(model, 10, 200)
    #pca_reach = ReachSet(pca_model)
    #stan_reach = ReachSet(stan_model)
    ross_reach = ReachSet(ross_pca)

    #pca_flow = mod_reach.computeReachSet(200, PCAStrat)
    #stan_flow = stan_reach.computeReachSet(200)
    ross_flow = ross_reach.computeReachSet(200, PCAStrat)
    sir_plot = Plot()
    #trajs = generate_traj(model, 10, 200)

    'Generate the trajectories and add them to the plot.'
    sir_plot.add(ross_flow)
    #sir_plot.add(pca_flow)
    sir_plot.plot(0,1,2)

    Timer.generate_stats()
