from models.quadcopter import Quadcopter, Quadcopter_PCA
from models.lotkavolterra import LotkaVolterra, LotkaVolterra_PCA
from models.rossler import Rossler, Rossler_PCA
from models.sir import SIR, SIR_PCA
from models.phos import Phosphorelay, Phosphorelay_PCA
from models.basic.basic import Basic
from models.basic.basic2 import Basic2

from kaa.reach import ReachSet
from kaa.timer import Timer
from settings import PlotSettings
from kaa.plotutil import Plot
from kaa.temp.pca_strat import PCAStrat

def test_all():

    num_steps = 100

    model_var_dict = {
        Basic : [0],
        Rossler: [0,1,2],
        SIR: [0,1,2],
        LotkaVolterra: [0,1,2],
        Quadcopter: [2,5,13],
        Phosphorelay: [0,1]
    }

    model_path_dict = {
        Basic : "/Users/edwardkim/Work/kaa-optimize/figures/Basic",
        Rossler: "/Users/edwardkim/Work/kaa-optimize/figures/Rossler",
        SIR: "/Users/edwardkim/Work/kaa-optimize/figures/SIR",
        LotkaVolterra: "/Users/edwardkim/Work/kaa-optimize/figures/LotVolt",
        Quadcopter: "/Users/edwardkim/Work/kaa-optimize/figures/Quad",
        Phosphorelay: "/Users/edwardkim/Work/kaa-optimize/figures/Phos"
    }

    for model, var_ind_list in model_var_dict.items():

        print("\n Generating Model: {}\n".format(str(model)))

        flowpipe = ReachSet(model()).computeReachSet(num_steps, PCAStrat)
        plot = Plot()
        plot.add(flowpipe)
        plot.plot(*var_ind_list, path=model_path_dict[model])



def test_all_pca():

    num_steps = 100

    model_var_dict = {
        Basic: [0],
        Rossler_PCA: [0,1,2],
        SIR_PCA: [0,1,2],
        LotkaVolterra_PCA: [0,1,2],
        Quadcopter_PCA: [2,5,13],
        Phosphorelay_PCA: [0,1]
    }

    model_path_dict = {
        Basic: "/Users/edwardkim/Work/kaa-optimize/figures/Basic",
        Rossler_PCA: "/Users/edwardkim/Work/kaa-optimize/figures/Rossler",
        SIR_PCA: "/Users/edwardkim/Work/kaa-optimize/figures/SIR",
        LotkaVolterra_PCA: "/Users/edwardkim/Work/kaa-optimize/figures/LotVolt",
        Quadcopter_PCA: "/Users/edwardkim/Work/kaa-optimize/figures/Quad",
        Phosphorelay_PCA: "/Users/edwardkim/Work/kaa-optimize/figures/Phos"
    }

    for model, var_ind_list in model_var_dict.items():

        print("\n Generating Model: {}\n".format(model))

        flowpipe = ReachSet(model()).computeReachSet(num_steps, PCAStrat)
        plot = Plot()
        plot.add(flowpipe)
        plot.plot(*var_ind_list, path=model_path_dict[model])
