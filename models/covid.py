from kaa.model import Model

class Covid(Model):


    def __init__(self, delta=0.1):

        sA, sI, A, I, Ra, Ri, D = sp.Symbol('sA'), sp.Symbol('sI'), sp.Symbol('A'), sp.Symbol('I'), sp.Symbol('Ra'), sp.Symbol('Ri'), sp.Symbol('D')

        dsA = sA + (-0.25 * sA * (A + I))*delta
        dsI = sI + (-0.25 * sI * (A + I))*delta
        dA = A + (0.25 * sA * (A + I) - gamma*A)*delta
        dI = I + (0.25 * sI * (A + I) - gamma*I)*delta
        dRa = Ra + (0.02*A)*delta
        dRi = Ri + (0.02*I)*delta
        dD = D + (0.02 * I)*delta

        vars = [sA, sI, A, I, Ra, Ri, D]
        dyns = [dsA, dsI, dA, dI, dRa, dRi, dD]
