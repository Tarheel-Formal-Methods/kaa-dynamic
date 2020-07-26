from kaa.pykodiak.pykodiak_interface import Kodiak

class Model:

    def __init__(self, bundle, f, vars):
        #Initial bundle
        self.bund = bundle
        self.f = f
        self.vars = vars

        for var in vars:
            Kodiak.add_variable(str(var))
