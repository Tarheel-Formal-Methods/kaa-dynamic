import sympy as sp
import numpy as np

from kaa.model import Model
from kaa.bundle import Bundle

'Lorentz Model'
class Lorentz_UnitBox(Model):

    def __init__(self, delta=1):

        x, y, z = sp.Symbol('x'), sp.Symbol('y')
        vars = [x, y]

        dim_sys = len(vars)

        dx = x + (10*(y-x))*delta
        dy = y + (x * (2.66 - z) - y)*delta
        dz = z + (x * y - 28*z)*delta

        dyns = [dx, dy, dz]

        num_direct = 3
        num_temps = 1

        L = np.zeros([num_direct, dim_sys])
        T = np.zeros([num_temps, dim_sys])

        L[0][0] = 1
        L[1][1] = 1
        L[2][2] = 1


        T[0][0] = 0; T[0][1] = 1; T[0][2] = 1;

        offu = np.zeros(num_direct)
        offl = np.zeros(num_direct)

        offu[0] = 15.001; offl[0] = -14.999;
        offu[1] = 15.001; offl[1] = -14.999;
        offu[2] = 36.001; offl[2] = -35.999;

        super().__init__(dyns, vars, T, L, offu, offl, name="Lorentz")
