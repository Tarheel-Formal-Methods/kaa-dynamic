import sympy as sp
import numpy as np

from kaa.model import Model
from kaa.reach import ReachSet
from kaa.plotutil import Plot
from settings import KaaSettings, PlotSettings

from kaa.timer import Timer

KaaSettings.SupressOutput = False
PlotSettings.save_fig = False

class HarOscRotate(Model):

    def __init__(self, delta=0.05):

        x,y = sp.Symbol('x'), sp.Symbol('y')

        dx = y
        dy = -x

        dyns  = [dx, dy]
        vars = [x, y]

        L = np.empty([2,2])
        T = np.empty([1,2])

        L[0] = [1, 1]
        L[1] = [-1, 1]

        T[0][0] = 0;
        T[0][1] = 1;

        #T[1][0] = 2
        #T[1][1] = 3

        offu = np.empty(2)
        offl = np.empty(2)

        offu[0] = 1
        offu[1] = 2
        #offu[2] = 5
        #offu[3] = 5

        offl[0] = 0
        offl[1] = -1.80
        #offl[2] = 5
        #offl[3] = 5
        #
        super().__init__(dyns, vars, T, L, offu, offl, name="HarOsc")


def test_haroscrotate():
    
    NUM_STEPS = 4

    model = HarOscRotate()
    mod_reach = ReachSet(model)

    mod_flow = mod_reach.computeReachSet(NUM_STEPS)

    SIR_PCA_ITER_STEPS = 1 #Number of steps between each recomputation of PCA Templates.
    'PCA Strategy Parameters'
    SIR_PCA_TRAJ_STEPS = 1 #Number of steps our sample trajectories should run.
    SIR_PCA_NUM_TRAJ = 100 #Number of sample trajectories we should use for the PCA routine.

    #pca_strat = PCAStrat(model, traj_steps=SIR_PCA_TRAJ_STEPS, num_trajs=SIR_PCA_NUM_TRAJ, iter_steps=SIR_PCA_ITER_STEPS)
    #mod_pca_flow = mod_reach.computeReachSet(NUM_STEPS, tempstrat=pca_strat)

    vdp_plot = Plot()
    vdp_plot.add(mod_flow, "HarOsc")
    #vdp_plot.add(mod_pca_flow, "HarOsc PCA")
    vdp_plot.plot2DPhase(0,1)

    Timer.generate_stats()
