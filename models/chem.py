import sympy as sp
import numpy as np

from kaa.model import Model
from kaa.bundle import Bundle

class Chem_UnitBox(Model):

    def __init__(self, beta, gamma, delta=0.2, init_box=((0.9999,1), (0.0000099, 0.00001), (0.0000099, 0.00001), (0.9999, 1))):

        x, y, z, s = sp.Symbol('x'), sp.Symbol('y'), sp.Symbol('z'), sp.Symbol('s')
        vars = (x, y, z, s)

        dim_sys = len(vars)

        dx = x + (-0.4*x + beta*y*z)*delta
        dy = y + (0.4*x - beta*y*z - gamma*y**2)*delta
        dz = z + (gamma*y**2)*delta
        ds = s + (x + y + z)*delta

        dyns = (dx, dy, dz, ds)

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

        super().__init__(dyns, vars, T, L, init_box, offl, offu, name="RobertsonChem")
