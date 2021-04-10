from sympy.utilities.lambdify import lambdify
from kaa.opts.kodiak import KodiakProd
from kaa.settings import KaaSettings
from kaa.bundle import Bundle

if KaaSettings.OptProd is KodiakProd:
    from kaa.pykodiak.pykodiak_interface import Kodiak

class Model:

    def __init__(self, f, vars, T, L, offu, offl, name="Model", compose=0):

        for _ in range(compose):
            var_sub = [ (var, f[var_idx]) for var_idx, var in enumerate(vars) ]
            f = [ func.subs(var_sub, simultaneous=True) for func in f ]

        'List of system dynamics.'
        self.f = f

        'List of system variables.'
        self.vars = vars

        'Dimension of system'
        self.dim = len(vars)

        'Name of system.'
        self.name = name

        'Initial bundle.'
        self.bund = Bundle(self, T, L, offu, offl)

        self.lambdified_f =  [lambdify(vars, f) for f in self.f]

    def __str__(self):
        return self.name
