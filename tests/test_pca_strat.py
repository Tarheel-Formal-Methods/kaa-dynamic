from kaa.reach import ReachSet
from kaa.plotutil import Plot
from models.sir import SIR, SIR_PCA

from kaa.timer import Timer
from kaa.experiutil import generate_traj
from kaa.temp.pca_strat import PCAStrat

def test_pca_strat():

    #pca_model = SIR_PCA()
    stan_model = SIR()
    
    #trajs = generate_traj(model, 10, 200)
    #pca_reach = ReachSet(pca_model)
    stan_reach = ReachSet(stan_model)

    #pca_flow = mod_reach.computeReachSet(200, PCAStrat)
    stan_flow = stan_reach.computeReachSet(200)

    sir_plot = Plot()
    #trajs = generate_traj(model, 10, 200)

    'Generate the trajectories and add them to the plot.'
    sir_plot.add(stan_flow)
    #sir_plot.add(pca_flow)
    sir_plot.plot(0,1,2)

    Timer.generate_stats()
