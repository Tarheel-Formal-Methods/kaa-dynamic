import numpy as np
import sympy as sp
import math

from kaa.model import Model

class HarOsc(Model):

    def __init__(self, delta=0.05):

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

        offu[0] = -4
        offu[1] = 1
        #offu[2] = 5
        #offu[3] = 5

        offl[0] = 5
        offl[1] = 0
        #offl[2] = 5
        #offl[3] = 5
        #
        super().__init__(dyns, vars, T, L, offu, offl, name="HarOsc")
