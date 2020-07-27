from kaa.pykodiak.pykodiak_interface import Kodiak

class Model:

    def __init__(self, bundle, f, vars, name="Model"):
        #Initial bundle
        self.bund = bundle
        self.f = f
        self.vars = vars
        self.name = name

        for var in vars:
            Kodiak.add_variable(str(var))
