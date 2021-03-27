import sympy as sp
import numpy as np

from kaa.model import Model
from kaa.bundle import Bundle

'Buckling Column Model'
class BuckCol_UnitBox(Model):

    def __init__(self, delta=0.05):
        x, y = sp.Symbol('x'), sp.Symbol('y')
        vars = [x, y]

        dim_sys = len(vars)

        dx = x + y*delta
        dy = y + (2*x - x**3 - 0.2*y + 0.1)*delta

        dyns = [dx, dy]

        num_direct = 2
        num_temps = 1

        L = np.zeros([num_direct, dim_sys])
        T = np.zeros([num_temps, dim_sys])

        L[0][0] = 1
        L[1][1] = 1

        T[0][0] = 0; T[0][1] = 1;

        offu = np.zeros(num_direct)
        offl = np.zeros(num_direct)

        offu[0] = -0.5; offl[0] = 0.4;
        offu[1] = -0.5; offl[1] = 0.4;

        super().__init__(dyns, vars, T, L, offu, offl, name="Buckling Column")
