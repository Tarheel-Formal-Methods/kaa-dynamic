import sympy as sp
import numpy as np

from kaa.model import Model
from kaa.bundle import Bundle

'CoupledVDP Model'
class CoupledVDP_UnitBox(Model):

    def __init__(self, delta=0.2):

        init_box = ((1.25, 1.55), (2.25, 2.35),(1.25, 1.55), (2.25, 2.35))
        x1, y1, x2, y2 = sp.Symbol('x1'), sp.Symbol('y1'), sp.Symbol('x2'), sp.Symbol('y2')
        vars = (x1, y1, x2, y2)

        dim_sys = len(vars)

        dx1 = x1 + (y1) * delta
        dy1 = y1 + ((1-x1**2)*y1 - 2*x1 + x2)*delta
        dx2 = x2 + (y2)*delta
        dy2 = y2 + ((1- x2**2)*y2 - x2 +(x1 - x2))*delta

        dyns = (dx1, dy1, dx2, dy2)

        num_direct = 4
        num_temps = 1

        L = np.zeros((num_direct, dim_sys))
        T = np.zeros((num_temps, dim_sys))

        L[0][0] = 1
        L[1][1] = 1
        L[2][2] = 1
        L[3][3] = 1

        T[0][0] = 0; T[0][1] = 1; T[0][2] = 2; T[0][3] = 3;

        offu = np.zeros(num_direct)
        offl = np.zeros(num_direct)

        offu[0] = 1.55; offl[0] = -1.25;
        offu[1] = 2.35; offl[1] = -2.25;
        offu[0] = 1.55; offl[0] = -1.25;
        offu[1] = 2.35; offl[1] = -2.25;

        super().__init__(dyns, vars, T, L, init_box, offl, offu, name="CoupledVDP")
