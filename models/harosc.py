import numpy as np
import sympy as sp
import math

from kaa.model import Model

class HarOsc(Model):

    def __init__(self, init_box=[[-5,-4],[0,1]]):

        x,y = sp.Symbol('x'), sp.Symbol('y')

        dx = 0.707107*x + 0.707107*y
        dy = -0.707107*x + 0.707107*y

        dyns  = [dx, dy]
        vars = [x, y]

        #L = np.empty([4,2])
        #T = np.empty([2,2])
        L = np.empty([2,2])
        T = np.empty([1,2])


        L[0] = [1, 0]
        L[1] = [0, 1]
        #L[0] = [1, 1]
        #L[1] = [-1, 1]

        T[0][0] = 0
        T[0][1] = 1

        #T[1][0] = 2
        #T[1][1] = 3

        offu = np.empty(2)
        offl = np.empty(2)

        for i in range(2):
            offu[i] = init_box[i][1]
            offl[i] = - init_box[i][0]


        super().__init__(dyns, vars, 0.1, T, L, init_box, offl, offu, name="HarOsc")
