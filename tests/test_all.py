from kaa.reach import ReachSet
from models.quadcopter import Quadcopter
from models.lotkavolterra import LotkaVolterra
from models.Rossler import Rossler
from models.SIR import SIR
from models.basic.basic import Basic
from models.basic.basic2 import Basic2

from kaa.timer import Timer
from kaa.settings import PlotSettings

PlotSettings.save_fig = True

def test_all():

    num_steps = 100

    model_var_dict = {
        Basic : [0],
        Rossler: [0,1,2],
        SIR: [0,1,2],
        LotkaVolterra: [0,1,2],
        Quadcopter: [2,5,13]
    }

    for model, var_ind_list in model_var_dict:

        flowpipe = ReachSet(model()).computeReachSet(num_steps)
        for var_ind in var_ind_list:
            flowpipe.plot2DProj(var_ind)
