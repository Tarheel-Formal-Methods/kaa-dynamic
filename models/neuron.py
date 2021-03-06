import sympy as sp
import numpy as np

from kaa.model import Model
from kaa.bundle import Bundle

'Fitz-Hugh-Nagumo Neuron Model'
class Neuron(Model):

    def __init__(self, delta=0.2, init_box=((0.9,1.1), (2.4,2.6))):

        x, y = sp.Symbol('x'), sp.Symbol('y')
        vars = (x, y)

        dim_sys = len(vars)

        dx = x + (x - x**3 - y + 0.875)*delta
        dy = y + (0.08*(x + 0.7 - 0.8*y))*delta

        dyns = (dx, dy)

        num_direct = 4
        num_temps = 6

        L = np.zeros((num_direct, dim_sys))
        T = np.zeros((num_temps, dim_sys))

        L[0][0] = 1
        L[1][1] = 1
        L[2][0] = -1; L[2][1] = 1;
        L[3][0] = 1; L[3][1] = 1;

        T[0][0] = 0; T[0][1] = 1;
        T[1][0] = 0; T[1][1] = 2;
        T[2][0] = 0; T[2][1] = 3;
        T[3][0] = 1; T[3][1] = 2;
        T[4][0] = 1; T[4][1] = 3;
        T[5][0] = 2; T[5][1] = 3; # SAme for neuron model, cvdp

        offu = np.zeros(num_direct)
        offl = np.zeros(num_direct)

        offu[2] = 10; offl[2] = 10;
        offu[3] = 10; offl[3] = 10;

        super().__init__(dyns, vars, T, L, init_box, offl, offu, name="Neuron")



'Fitz-Hugh-Nagumo Neuron Model'
class Neuron_UnitBox(Model):

    def __init__(self, delta=0.2, init_box=((0.9,1.1), (2.4,2.6))):

        x, y = sp.Symbol('x'), sp.Symbol('y')
        vars = (x, y)

        dim_sys = len(vars)

        dx = x + (x - x**3 - y + 0.875)*delta
        dy = y + (0.08*(x + 0.7 - 0.8*y))*delta

        dyns = (dx, dy)

        num_direct = 2
        num_temps = 1

        L = np.zeros((num_direct, dim_sys))
        T = np.zeros((num_temps, dim_sys))

        L[0][0] = 1
        L[1][1] = 1

        T[0][0] = 0; T[0][1] = 1;

        offu = np.zeros(num_direct)
        offl = np.zeros(num_direct)

        super().__init__(dyns, vars, T, L, init_box, offl, offu, name="Neuron")
