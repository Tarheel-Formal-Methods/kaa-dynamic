import sympy as sp
import numpy as np

from kaa.model import Model
from kaa.bundle import Bundle

class CoupledVDP(Model):

    def __init__(self, delta=0.1, init_box=((1.25,1.55), (2.25,2.35),(1.25,1.55), (2.25,2.35)), mu=1):

        x1, y1, x2, y2  = sp.Symbol('x1'), sp.Symbol('y1'), sp.Symbol('x2'), sp.Symbol('y2')
        vars = (x1, y1, x2, y2)

        dim_sys = len(vars)

        dx1 = x1 + y1*delta
        dy1 = y1 + (mu*(1-x1**2)*y1 - 2*x1 + x2)*delta
        dx2 = x2 + y2*delta
        dy2 = y2 + (mu*(1-x2**2)*y2 - 2*x2 + x1)*delta

        dyns = (dx1, dy1, dx2, dy2)

        num_direct = 9
        num_temps = 4

        L = np.zeros((num_direct, dim_sys))
        T = np.zeros((num_temps, dim_sys))

        for i in range(4):
            L[i][i] = 1

        L[4][0] = 1; L[4][1] = 1;
        L[5][2] = 1; L[5][3] = 1;
        L[6][0] = 1; L[6][2] = 1;
        L[7][1] = 1; L[7][3] = 1;
        L[8][0] = 1; L[8][1] = 1; L[8][2] = -1;

        T[0][0] = 0; T[0][1] = 1; T[0][2] = 2; T[0][3] = 3;
        T[1][0] = 1; T[1][1] = 2; T[1][2] = 3; T[1][3] = 4;
        T[2][0] = 4; T[2][1] = 5; T[2][2] = 6; T[2][3] = 8;
        T[3][0] = 3; T[3][1] = 4; T[3][2] = 5; T[3][3] = 6;
        #T[4][0] = 4; T[4][1] = 5; T[4][2] = 6; T[4][3] = 7;
        #T[5][0] = 4; T[5][1] = 5; T[5][2] = 7; T[5][3] = 8;

        offu = np.zeros(num_direct)
        offl = np.zeros(num_direct)

        offu[4] = 3; offl[4] = 3;
        offu[5] = 3; offl[5] = 3;
        offu[6] = 3; offl[6] = 3;
        offu[7] = 3; offl[7] = 3;
        offu[8] = 3; offl[8] = 3;

        super().__init__(dyns, vars, delta, T, L, init_box, offl, offu, name="CoupledVDP")


class CoupledVDP_UnitBox(Model):

    def __init__(self, delta=0.2, init_box=((1.25,1.55), (2.25,2.35),(1.25,1.55), (2.25,2.35), (60,80)), mu=1):

        x1, y1, x2, y2, b = sp.Symbol('x1'), sp.Symbol('y1'), sp.Symbol('x2'), sp.Symbol('y2'), sp.Symbol('b')
        vars = (x1, y1, x2, y2, b)

        dim_sys = len(vars)

        dx1 = x1 + y1*delta
        dy1 = y1 + (mu*(1-x1**2)*y1 - b*(x2 - x1) - x1)*delta
        dx2 = x2 + y2*delta
        dy2 = y2 + (mu*(1-x2**2)*y2 - b*(x2 - x1) - x2)*delta
        db = b


        dyns = (dx1, dy1, dx2, dy2, db)

        num_direct = 5
        num_temps = 1

        L = np.zeros((num_direct, dim_sys))
        T = np.zeros((num_temps, dim_sys))

        L[0][0] = 1
        L[1][1] = 1
        L[2][2] = 1
        L[3][3] = 1
        L[4][4] = 1

        T[0][0] = 0; T[0][1] = 1; T[0][2] = 2; T[0][3] = 3; T[0][4] = 4;

        offu = np.zeros(num_direct)
        offl = np.zeros(num_direct)

        super().__init__(dyns, vars, delta, T, L, init_box, offl, offu, name="CoupledVDP")
