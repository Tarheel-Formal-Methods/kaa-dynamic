'''
Python interface file for pykodiak.so

Staley Bak
Oct 2018
'''

import ctypes
import os

import numpy as np
from numpy.ctypeslib import ndpointer

from sympy import Mul, Expr, Add, Pow, Symbol, Number, sin, cos, atan
from sympy.parsing.sympy_parser import parse_expr

def get_script_path(filename):
    '''get the path this script, pass in __file__ for the filename'''
    return os.path.dirname(os.path.realpath(filename))

class Kodiak:

    def __init__(self):
        self.initialized = False
        self.variables = {}
        self.init()

    def init(self):
        'initialize the static members'

        if not self.initialized:
            self.initialized = True

            lib_path = os.path.join(get_script_path(__file__), 'pykodiak.so')
            self.lib = ctypes.CDLL(lib_path)

            self._init = self.lib.init
            self._init.restype = None
            self._init.argtypes = []

            self._use_bernstein = self.lib.useBernstein
            self._use_bernstein.restype = None
            self._use_bernstein.argtypes = [ctypes.c_int]

            self._set_precision = self.lib.setPrecision
            self._set_precision.restype = None
            self._set_precision.argtypes = [ctypes.c_int]

            self._add_variable = self.lib.addVariable
            self._add_variable.restype = None
            self._add_variable.argtypes = [ctypes.c_char_p]

            self._lookup_variable = self.lib.lookupVariable
            self._lookup_variable.restype = ctypes.c_int
            self._lookup_variable.argtypes = [ctypes.c_char_p]

            self._print_expression = self.lib.printExpression
            self._print_expression.restype = None
            self._print_expression.argtypes = [ctypes.c_int]

            self._make_double = self.lib.makeDouble
            self._make_double.restype = ctypes.c_int
            self._make_double.argtypes = [ctypes.c_double]

            self._make_mult = self.lib.makeMult
            self._make_mult.restype = ctypes.c_int
            self._make_mult.argtypes = [ctypes.c_int, ctypes.c_int]

            self._make_add = self.lib.makeAdd
            self._make_add.restype = ctypes.c_int
            self._make_add.argtypes = [ctypes.c_int, ctypes.c_int]

            self._make_sq = self.lib.makeSq
            self._make_sq.restype = ctypes.c_int
            self._make_sq.argtypes = [ctypes.c_int]

            self._make_sqrt = self.lib.makeSqrt
            self._make_sqrt.restype = ctypes.c_int
            self._make_sqrt.argtypes = [ctypes.c_int]

            self._make_intpow = self.lib.makeIntPow
            self._make_intpow.restype = ctypes.c_int
            self._make_intpow.argtypes = [ctypes.c_int, ctypes.c_int]

            self._make_sin = self.lib.makeSin
            self._make_sin.restype = ctypes.c_int
            self._make_sin.argtypes = [ctypes.c_int]

            self._make_cos = self.lib.makeCos
            self._make_cos.restype = ctypes.c_int
            self._make_cos.argtypes = [ctypes.c_int]

            self._make_atan = self.lib.makeAtan
            self._make_atan.restype = ctypes.c_int
            self._make_atan.argtypes = [ctypes.c_int]

            self._free_stack = self.lib.free_stack
            self._free_stack.restype = None
            self._free_stack.argtypes = []

            self._minmax_diff = self.lib.minmax_diff
            self._minmax_diff.restype = None
            self._minmax_diff.argtypes = [ctypes.c_int, # nonlinear expression
                                         ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), ctypes.c_int,
                                         ctypes.c_double, # bias
                                         ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), ctypes.c_int,
                                         ctypes.c_int,
                                         ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), ctypes.c_int]



            #### initialize with default values ####
            self._init()
            self.reset()

    def reset(self):
        '''reset settings to default'''

        self.use_bernstein(True)
        self.set_precision(-6) # 10^-6 precision

    def to_expression(self, str_exp):
        'convert a string expression to an Kodiak expression object'

        return self.sympy_to_kodiak(parse_expr(str_exp))


    def minmax_diff(self, nonlinear_exp_index, linear_approx, linear_bias, bounds):
        ''''return the lower and upper bound of the passed-in nonlinear function minus the linear approximation,
        within the passed-in bounds'''

        if not isinstance(linear_approx, np.ndarray):
            linear_approx = np.array(linear_approx, dtype=float)

        if not isinstance(bounds, np.ndarray):
            bounds = np.array(bounds, dtype=float)

        assert len(linear_approx) == len(bounds), f"Length of LinApp: {len(bounds)}, Length of Bounds: {len(linear_approx)}"

        # rv = np.array([0, 0], dtype=float)

        rv = np.array([0, 0, 0, 0], dtype=float)

        self._minmax_diff(nonlinear_exp_index, linear_approx, len(linear_approx), linear_bias,
                         bounds, bounds.shape[0], bounds.shape[1], rv, rv.shape[0])

        #todo change back
        #return rv[0], rv[1]
        return rv[0], rv[1], rv[2], rv[3]

    def free_stack(self):
        self._free_stack()

    def use_bernstein(self, use_bernstein):
        'should we use bernstein polynomials for optimization (false = interval arithmetic), default: True'
        self._use_bernstein(1 if use_bernstein else 0)


    def set_precision(self, precision):
        'set the prevision for the answer, more negative = more accurage, default: -9 (10^-9)'
        self._set_precision(precision)


    def add_variable(self, name):
        'add an optimization variable (should be called in order)'

        #self.init()
        self.variables[name] = len(self.variables)

        n = name.encode('utf-8')
        self._add_variable(n)


    def lookup_variable(self, name):
        'lookup the expression index for a variable previously-added with add_variable()'

        #self.init()

        assert name in self.variables, f"variable '{name}' must be added first with Kodiak.add_variable()"

        n = name.encode('utf-8')
        return self._lookup_variable(n)


    def print_expression(self, i):
        'print the expression with the given index to stdout'

        #self.init()
        self._print_expression(i)


    def make_double(self, d):
        'make and return an expression index for a double number'

        #self.init()
        return self._make_double(d)


    def make_mult(self, a, b):
        'make and return an expression index for a multiplication of the two passed-in expression indices'

        #self.init()
        return self._make_mult(a, b)


    def make_add(self, a, b):
        'make and return an expression index for an addition of the two passed-in expression indices'

        #self.init()
        return self._make_add(a, b)


    def make_sq(self, a):
        'make and return an expression index for the square of the the passed-in expression index'

        #self.init()
        return self._make_sq(a)


    def make_sqrt(self, a):
        'make and return an expression index for the square of the the passed-in expression index'

        #self.init()
        return self._make_sqrt(a)


    def make_intpow(self, a, intb):
        '''make and return an expression index for passed-in expression index raised to the integer intb
        note: intb is an integer, not an expression index
        '''

        #self.init()
        return self._make_intpow(a, intb)


    def make_sin(self, a):
        'make and return an expression index for the sine of the the passed-in expression index'

        #self.init()
        return self._make_sin(a)


    def make_cos(self, a):
        'make and return an expression index for the cosine of the the passed-in expression index'

        #self.init()
        return self._make_cos(a)


    def make_atan(self, a):
        'make and return an expression index for the atan of the the passed-in expression index'

        #self.init()
        return self._make_atan(a)


    def sympy_to_kodiak(self, sympy_exp):
        '''convert a sympy expression to Kodiak expression

        this function actually returns an int, which is the expression index in the c++ 'reals' vector that
        represents the (newly-constructed) expression
        '''

        #self.init()

        rv = None
        e = sympy_exp

        if not isinstance(e, Expr):
            raise RuntimeError("Expected sympy Expr: " + repr(e))

        if isinstance(e, Symbol):
            rv = self.lookup_variable(e.name)

            if rv is None:
                raise RuntimeError("No var was corresponds to symbol '" + str(e) + "'")
        elif isinstance(e, Number):
            rv = self.make_double(float(e))
        elif isinstance(e, Mul):
            rv = self.sympy_to_kodiak(e.args[0])

            for arg in e.args[1:]:
                val = self.sympy_to_kodiak(arg)
                rv = self.make_mult(rv, val)
        elif isinstance(e, Add):
            rv = self.sympy_to_kodiak(e.args[0])

            for arg in e.args[1:]:
                val = self.sympy_to_kodiak(arg)
                rv = self.make_add(rv, val)
        elif isinstance(e, Pow):
            term = self.sympy_to_kodiak(e.args[0])
            exponent = e.args[1]

            assert isinstance(exponent, Number), f"exponent must be a number: {e}, (type {type(e)})"

            if float(exponent) == 2:
                # squared
                rv = self.make_sq(term)
            elif float(exponent) == 0.5:
                # sqrt
                rv = self.make_sqrt(term)
            else:
                assert float(exponent) == int(exponent), f"exponent must be an integer (or 0.5): {e}"

                # intpow
                rv = self.make_intpow(term, exponent)

            # non-integer powers are not supported in Kodiak
        elif isinstance(e, sin):
            a = self.sympy_to_kodiak(e.args[0])
            rv = self.make_sin(a)
        elif isinstance(e, cos):
            a = self.sympy_to_kodiak(e.args[0])
            rv = self.make_cos(a)
        elif isinstance(e, atan):
            a = self.sympy_to_kodiak(e.args[0])
            rv = self.make_atan(a)

        # for all expression supported by Kodiak, see Real.hpp

        assert rv is not None, f"conversion of '{e}' to a Kodiak expression unsupported (type {type(e)})"

        return rv
