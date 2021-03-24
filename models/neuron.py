import sympy as sp
import numpy as np

from kaa.model import Model
from kaa.bundle import Bundle

'Fitz-Hugh-Nagumo Neuron Model'
class Neuron_UnitBox(Model):

    def __init__(self, delta=0.2):

        x, y = sp.Symbol('x'), sp.Symbol('y')
        vars = [x, y]

        dim_sys = len(vars)

        dx1 = x + (x - x**3 - y + 0.875)*delta
        dx2 = y + (0.08*(x + 0.7 - 0.8*y))*delta

        dyns = [x, y]

        num_direct = 2
        num_temps = 1

        L = np.zeros([num_direct, dim_sys])
        T = np.zeros([num_temps, dim_sys])

        L[0][0] = 1
        L[1][1] = 1

        T[0][0] = 0; T[0][1] = 1;

        offu = np.zeros(num_direct)
        offl = np.zeros(num_direct)

        offu[0] = 1.1; offl[0] = -0.9;
        offu[1] = 2.6; offl[1] = -2.4;

        super().__init__(dyns, vars, T, L, offu, offl, name="Neuron")
