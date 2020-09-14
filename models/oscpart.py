import numpy as np
import sympy as sp

from kaa.model import Model

class OscPart(Model):

    def __init__(self):

        x, y, z = sp.Symbol('x'), sp.Symbol('y'), sp.Symbol('z')

        vars = [x, y, z]

        dx = -0.1 * x - y
        dy = x - 0.1 * y
        dz = -0.15 * z

        dyns  = [dx, dy, dz]

        L = np.empty([3,3])
        T = np.empty([1,3])

        L[0] = [1, 0, 0]
        L[1] = [0, 1, 0]
        L[2] = [0, 0, 1]

        T[0][0] = 0
        T[0][1] = 1
        T[0][2] = 2

        offu = np.empty(3)
        offl = np.empty(3)

        offu[0] = 0.1
        offu[1] = 1
        offu[2] = 1

        offl[0] = 0.1
        offl[1] = -0.8
        offl[2] = -0.9

        super().__init__(dyns, vars, T, L, offu, offl, name="OscPart")
