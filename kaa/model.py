from sympy.utilities.lambdify import lambdify
from kaa.opts.kodiak import KodiakProd
from settings import KaaSettings
from kaa.bundle import Bundle

if KaaSettings.OptProd is KodiakProd:
    pass

class Model:

    def __init__(self, f, vars, step_size, T, L, init_box, offl, offu, name="Model", compose=0):

        for _ in range(compose):
            var_sub = [(var, f[var_idx]) for var_idx, var in enumerate(vars)]
            f = [func.subs(var_sub, simultaneous=True) for func in f]

        'List of system dynamics.'
        self.f = f

        'List of system variables.'
        self.vars = vars

        'Step size'
        self.step_size = step_size

        'Dimension of system'
        self.dim = len(vars)

        'Name of system.'
        self.name = name

        self.__set_init_box(init_box, offl, offu)

        self.init_box = init_box

        'Initial bundle.'
        self.bund = Bundle(self, T, L, offu, offl)

        if not KaaSettings.Parallelize:
            self.lambdified_f = lambdify(vars, self.f)

    def __str__(self):
        return self.name

    def __set_init_box(self, init_box, offl, offu):
        assert len(init_box) == self.dim, f"dimensions of init box have to same as that of system's, Init Box Dim: {len(init_box)}, Sys Dim: {self.dim}"

        for idx in range(self.dim):
            offl[idx] = -init_box[idx][0]
            offu[idx] = init_box[idx][1]
