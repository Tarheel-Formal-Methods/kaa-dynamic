import numpy as np
import sympy as sp

from kaa.model import Model

class InvertPend(Model):

    def __init__(self, delta = 0.1, init_box=((0.25,0.3),(0.25, 0.3))):

        x,y = sp.Symbol('x'), sp.Symbol('y')

        dx = x + y * delta
        dy = y + (-0.2*y - sp.sin(x)) * delta

        dyns  = [dx, dy]
        vars = [x, y]

        dim_sys = 2
        num_dirs = 4
        num_temps = 6

        L = np.zeros((num_dirs, dim_sys))

        L[0][0] = 1
        L[1][1] = 1
        L[2][0] = -1; L[2][1] = 1;
        L[3][0] = 1; L[3][1] = 1;

        T = np.zeros((num_temps, dim_sys))
        T[0][0] = 0; T[0][1] = 1;
        T[1][0] = 0; T[1][1] = 2;
        T[2][0] = 0; T[2][1] = 3;
        T[3][0] = 1; T[3][1] = 2;
        T[4][0] = 1; T[4][1] = 3;
        T[5][0] = 2; T[5][1] = 3;

        offu = np.zeros(num_dirs);
        offl = np.zeros(num_dirs);

        offu[2] = 10; offl[2] = 10;
        offu[3] = 10; offl[3] = 10;

        super().__init__(dyns, vars, delta, T, L, init_box, offl, offu, name="InvertPend")


class InvertPend_UnitBox(Model):

    def __init__(self, delta = 0.1, init_box=((0.25,0.3),(0.25, 0.3))):

        x,y = sp.Symbol('x'), sp.Symbol('y')

        dx = x + y * delta
        dy = y + (-0.2*y - sp.sin(x)) * delta
        ((0.25, 0.3), (0.25, 0.3))
        vars = (x,y)
        dyns = (dx, dy)

        dim_sys = 2
        num_dirs = 2
        num_temps = 1

        L = np.zeros((num_dirs, dim_sys))

        L[0][0] = 1
        L[1][1] = 1

        T = np.zeros((num_temps, dim_sys))
        T[0][0] = 0; T[0][1] = 1;

        offu = np.zeros(num_dirs);
        offl = np.zeros(num_dirs);

        super().__init__(dyns, vars, delta, T, L, init_box, offl, offu, name="InvertPend")
