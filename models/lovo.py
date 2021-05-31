import sympy as sp
import numpy as np

from kaa.model import Model
from kaa.bundle import Bundle

class LOVO_UnitBox(Model):

    def __init__(self, delta=0.05, init_box=((1.288, 1.312), (0.99999, 1))):

        x, y = sp.Symbol('x'), sp.Symbol('y')
        vars = (x, y)

        dim_sys = len(vars)

        dx = x + (3*x - 3*x*y)*delta
        dy = y + (x*y - y)*delta

        dyns = (dx, dy)

        num_direct = 2
        num_temps = 1

        L = np.zeros((num_direct, dim_sys))
        T = np.zeros((num_temps, dim_sys))

        L[0][0] = 1
        L[1][1] = 1

        T[0][0] = 0; T[0][1] = 1;

        offu = np.zeros(num_direct)
        offl = np.zeros(num_direct)

        super().__init__(dyns, vars, T, L, init_box, offl, offu, name="LOVO21")
