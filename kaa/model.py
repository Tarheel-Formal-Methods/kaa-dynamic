from kaa.opts.kodiak import KodiakProd
from kaa.settings import KaaSettings
from kaa.bundle import Bundle

if KaaSettings.OptProd is KodiakProd:
    from kaa.pykodiak.pykodiak_interface import Kodiak

class Model:

    def __init__(self, f, vars, T, L, offu, offl, name="Model"):
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

        if KaaSettings.OptProd is KodiakProd:
            for var in self.vars:
                Kodiak.add_variable(str(var))

    def __str__(self):
        return self.name
