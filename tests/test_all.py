from models.quadcopter import Quadcopter
from models.lotkavolterra import LotkaVolterra
from models.rossler import Rossler
from models.sir import SIR
from models.phos import Phosphorelay
from models.basic.basic import Basic
from models.basic.basic2 import Basic2

from kaa.reach import ReachSet
from kaa.timer import Timer
from kaa.settings import PlotSettings
from kaa.plotutil import Plot

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

        print("\n Generating Model: {}\n".format(model))

        flowpipe = ReachSet(model()).computeReachSet(num_steps)
        plot = Plot()
        plot.add(flowpipe)
        plot.plot(*var_ind_list, path=model_path_dict[model])
