import sympy as sp
import numpy as np

from kaa.model import Model
from kaa.bundle import Bundle

class Duffing_UnitBox(Model):

    def __init__(self, delta=0.05):

        x1, x2 = sp.Symbol('x1'), sp.Symbol('x2')
        vars = [x1, x2]

        dim_sys = len(vars)

        dx1 = x1 + x2*delta
        dx2 = x2 + (-x1 - 2*0.3*x_2 - x1**3)*delta

        dyns = [dx1, dx2]

        num_direct = 2
        num_temps = 1

        L = np.zeros([num_direct, dim_sys])
        T = np.zeros([num_temps, dim_sys])

        L[0][0] = 1
        L[1][1] = 1

        T[0][0] = 0; T[0][1] = 1;

        offu = np.zeros(num_direct)
        offl = np.zeros(num_direct)

        offu[0] = 2.51; offl[0] = -2.49;
        offu[1] = 1.51; offl[1] = -1.49;

        super().__init__(dyns, vars, T, L, offu, offl, name="Duffing")
