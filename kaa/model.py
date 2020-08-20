from kaa.opts.kodiak import KodiakProd
from kaa.settings import KaaSettings

if isinstance(KaaSettings.OptProd, KodiakProd):
    from kaa.pykodiak.pykodiak_interface import Kodiak

class Model:

    def __init__(self, bundle, f, vars, name="Model"):
        'Initial bundle.'
        self.bund = bundle

        'List of system dynamics.'
        self.f = f

        'List of system variables.'
        self.vars = vars

        'Name of system.'
        self.name = name

        if isinstance(KaaSettings.OptProd, KodiakProd):
            for var in self.vars:
                Kodiak.add_variable(str(var))
