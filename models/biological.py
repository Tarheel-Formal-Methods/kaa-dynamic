class Biological_UnitBox(Model):

    def __init__(self, delta=0.1):

        init_box = ((0.99, 1.01), (0.99, 1.01), (0.99, 1.01), (0.99, 1.01), (0.99, 1.01), (0.99, 1.01), (0.99, 1.01))

        x1, x2, x3, x4, x5, x6, x7 = sp.Symbol('x1'), sp.Symbol('x2'), sp.Symbol('x3'), sp.Symbol('x4'), sp.Symbol('x5'), sp.Symbol('x6'), sp.Symbol('x7')

        dx1 = x1 + (-0.4 * x1 + 5*x3*x4)*delta
        dx2 = x2 + (0.4*x1 - x2)*delta
        dx3 = x3 + (x2 - 5*x3*x4)*delta
        dx4 = x4 + (5*x5*x6 - 5*x3*x4)*delta
        dx5 = x5 + (-5*x5*x6 + 5*x3*x4)*delta
        dx6 = x6 + (0.5*x7 - 5*x5*x6)*delta
        dx7 = x7 + (-0.5*x7 + 5*x5*x6)*delta

        vars = (x1, x2, x3, x4, x5, x6, x7)
        dyns = (dx1, dx2, dx3, dx4, dx5, dx6, dx7)
        sys_dim = len(vars)

        num_dirs = 7
        num_temps = 1


        L = np.zeros((num_dirs, sys_dim))
        for i in range(sys_dim):
            L[i][i] = 1

        T = np.zeros((num_temps, sys_dim))
        T[0][0] = 0; T[0][1] = 1; T[0][2] = 2; T[0][3] = 3; T[0][4] = 4; T[0][5] = 5; T[0][6] = 6;

        offu = np.zeros(num_dirs)
        offl = np.zeros(num_dirs)

        super().__init__(dyns, vars, T, L, init_box, offl, offu, name="Biological")
