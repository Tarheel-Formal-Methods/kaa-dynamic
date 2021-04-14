import sympy as sp
import numpy as np

from kaa.model import Model
from kaa.bundle import Bundle

'JetEngine Model'
class JetEngine_UnitBox(Model):

    def __init__(self, delta=0.2, init_box=((0.8,1.2), (0.8,1.2))):

        x, y = sp.Symbol('x'), sp.Symbol('y')
        vars = [x, y]

        dim_sys = len(vars)

        dx = x + (-y - 1.5*x**2 - 0.5*x**3 - 0.5)*delta
        dy = y + (3*x - y)*delta

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

        super().__init__(dyns, vars, T, L, init_box, offl, offu, name="JetEngine")
