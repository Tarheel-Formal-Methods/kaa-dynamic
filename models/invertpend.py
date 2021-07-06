import numpy as np
import sympy as sp

from kaa.model import Model

class InvertPend(Model):

    def __init__(self, delta = 0.1, init_box=[[0,1],[0,1]]):

        x,y = sp.Symbol('x'), sp.Symbol('y')

        dx = x + y * delta
        dy = y + (-0.2*y - sp.sin(x)) * delta

        dyns  = [dx, dy]
        vars = [x, y]

        L = np.empty((2,2))
        T = np.empty((1,2))


        L[0] = [1, 0]
        L[1] = [0, 1]

        T[0][0] = 0
        T[0][1] = 1

        offu = np.empty(2)
        offl = np.empty(2)

        for i in range(2):
            offu[i] = init_box[i][1]
            offl[i] = - init_box[i][0]


        super().__init__(dyns, vars, delta, T, L, init_box, offl, offu, name="HarOsc")
