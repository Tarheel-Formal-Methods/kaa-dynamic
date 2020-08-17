from kaa.pykodiak.pykodiak_interface import Kodiak

class Model:

    def __init__(self, bundle, f, vars, name="Model", initial_set=None):
        'Initial bundle.'
        self.bund = bundle

        'List of system dynamics.'
        self.f = f

        'List of system variables.'
        self.vars = vars

        'Name of system.'
        self.name = name

        'Intevals reprsenting initial set in the order of self.vars above.'
        self.initial_set = initial_set

        'Ensuring Kodiak works as it should for us.'
        for var in vars:
            Kodiak.add_variable(str(var))
