import sympy as sp
import numpy as np

from kaa.pykodiak.pykodiak_interface import Kodiak

def test_seg_fault():
    x, y = sp.Symbol('x'), sp.Symbol('y')
    poly = -134960909.098539*x + 82082638596.6177*y - 4.65914220457606e-19*(1 - 6.81812756737304e+21*(-0.000785800134345946*x + y - 0.0683305085482289)**2)*(-875948208.091513*x - 6116441028.06391*y + 500352654.2844) - 5602155388.02084
    kodiak = Kodiak()

    kodiak.add_variable(str(x))
    kodiak.add_variable(str(y))
    kodiak_poly = kodiak.sympy_to_kodiak(poly)


    bounds = [[0,1] for _ in range(2)]
    jac_mat = np.zeros((2, 2))

    #print("Calling Kodiak")
    #print(f"INPUT POLY: {self.poly}")
    lb, ub, _, _ = kodiak.minmax_diff(kodiak_poly, jac_mat, 0, bounds)
    #print("Out of Kodiak")
