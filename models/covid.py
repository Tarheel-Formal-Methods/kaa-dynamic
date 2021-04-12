import sympy as sp
import numpy as np

from kaa.model import Model


class Covid(Model):

    def __init__(self, delta=0.1):

        sA, sI, A, I, Ra, Ri, D = sp.Symbol('sA'), sp.Symbol('sI'), sp.Symbol('A'), sp.Symbol('I'), sp.Symbol('Ra'), sp.Symbol('Ri'), sp.Symbol('D')

        dsA = sA + (-0.25 * sA * (A + I))*delta
        dsI = sI + (-0.25 * sI * (A + I))*delta
        dA = A + (0.25 * sA * (A + I) - 0.02*A)*delta
        dI = I + (0.25 * sI * (A + I) - 0.02*I)*delta
        dRa = Ra + (0.02*A)*delta
        dRi = Ri + (0.02*I)*delta
        dD = D + (0.02 * I)*delta

        vars = [sA, sI, A, I, Ra, Ri, D]
        dyns = [dsA, dsI, dA, dI, dRa, dRi, dD]
        sys_dim = len(vars)

        num_dirs = 10
        num_temps = 4

        L = np.zeros((num_dirs, sys_dim))
        for i in range(sys_dim):
            L[i][i] = 1

        L[7][2] = 1; L[7][3] = 1;
        L[8][4] = 1; L[8][5] = 1;
        L[9][2] = 1; L[9][3] = 1; L[9][4] = 1; L[9][5] = 1;

        T = np.zeros([num_temps,sys_dim])
        T[0][0] = 0; T[0][1] = 1; T[0][2] = 2; T[0][3] = 3; T[0][4] = 4; T[0][5] = 5; T[0][6] = 6;
        T[1][0] = 0; T[1][1] = 1; T[1][2] = 2; T[1][3] = 7; T[1][4] = 4; T[1][5] = 5; T[1][6] = 6;
        T[2][0] = 0; T[2][1] = 1; T[2][2] = 2; T[2][3] = 7; T[2][4] = 8; T[2][5] = 5; T[2][6] = 6;
        T[3][0] = 0; T[3][1] = 1; T[3][2] = 2; T[3][3] = 7; T[3][4] = 9; T[3][5] = 5; T[3][6] = 6;

        offu = np.zeros(num_dirs)
        offl = np.zeros(num_dirs)

        offu[0] = 0.7; offl[0] = -0.69;
        offu[1] = 0.1; offl[1] = -0.09;
        offu[2] = 0.15; offl[2] = -0.14;
        offu[3] = 0.05; offl[3] = -0.04;
        offu[4] = 0.001; offl[4] = -0.00099;
        offu[5] = 0.001; offl[5] = -0.00099;
        offu[6] = 0.001; offl[6] = -0.00099;

        super().__init__(dyns, vars, T, L, offu, offl, name="covid")

class Covid_UnitBox(Model):

    def __init__(self, delta=0.1):

        sA, sI, A, I, Ra, Ri, D = sp.Symbol('sA'), sp.Symbol('sI'), sp.Symbol('A'), sp.Symbol('I'), sp.Symbol('Ra'), sp.Symbol('Ri'), sp.Symbol('D')

        dsA = sA + (-0.25 * sA * (A + I))*delta
        dsI = sI + (-0.25 * sI * (A + I))*delta
        dA = A + (0.25 * sA * (A + I) - 0.02*A)*delta
        dI = I + (0.25 * sI * (A + I) - 0.02*I)*delta
        dRa = Ra + (0.02*A)*delta
        dRi = Ri + (0.02*I)*delta
        dD = D + (0.02 * I)*delta

        vars = [sA, sI, A, I, Ra, Ri, D]
        dyns = [dsA, dsI, dA, dI, dRa, dRi, dD]
        sys_dim = len(vars)

        num_dirs = 7
        num_temps = 1


        L = np.zeros([num_dirs, sys_dim])
        for i in range(sys_dim):
            L[i][i] = 1

        T = np.zeros([num_temps,sys_dim])
        T[0][0] = 0; T[0][1] = 1; T[0][2] = 2; T[0][3] = 3; T[0][4] = 4; T[0][5] = 5; T[0][6] = 6;

        offu = np.zeros(num_dirs)
        offl = np.zeros(num_dirs)

        offu[0] = 0.7; offl[0] = -0.69;
        offu[1] = 0.1; offl[1] = -0.09;
        offu[2] = 0.15; offl[2] = -0.14;
        offu[3] = 0.05; offl[3] = -0.04;
        offu[4] = 0.001; offl[4] = -0.00099;
        offu[5] = 0.001; offl[5] = -0.00099;
        offu[6] = 0.001; offl[6] = -0.00099;

        super().__init__(dyns, vars, T, L, offu, offl, name="covid")
