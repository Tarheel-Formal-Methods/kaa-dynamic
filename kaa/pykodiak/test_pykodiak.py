'''
Test / demo code for pykodiak interface.

Stanley Bak
Oct 2018
'''

from pykodiak_interface import Kodiak
from sympy.parsing.sympy_parser import parse_expr

def setup_function(_):
    '''setup method for each test'''

    #Kodiak.use_bernstein(True)
    #Kodiak.set_precision(-3) # 10^-3 precision

def test1():
    'main demo code'

    derivatives = ["x1", "(1-x0**2)*x1-x0"]

     # create sympy expressions from strings
    sympy_ders = [parse_expr(d) for d in derivatives]

     # add variables, in order, to kodiak
    Kodiak.add_variable('x0')
    Kodiak.add_variable('x1')

     # convert sympy expressions into kodiak expressions
    kodiak_ders = [Kodiak.sympy_to_kodiak(d) for d in sympy_ders]

    jac_mat = [[0., 1.], [-7.44, -0.96]]
    bias = 9.016
    bounds = [[1.25, 1.55], [2.28, 2.32]]

    start = time.time()

    iterations = 1000
    for _ in range(iterations):
        lb, ub, _, _ = Kodiak.minmax(kodiak_ders[1], jac_mat[1], bias, bounds)

    print("runtime for {} iterations: {} sec".format(iterations, round(time.time() - start, 3)))
    print("Kodiak computed enclosure: [{}, {}]".format(lb, ub))

def test_enclosure():
    Kodiak.add_variable('x0')
    Kodiak.add_variable('x1')
    Kodiak.add_variable('x2')
    
    expression = Kodiak.to_expression("(1-x0**2)*x1-x0 + (x2 - x0)")

    jac_mat = [-8.44, -0.96, 1.]
    bias = 9.016
    bounds = [[1.24999, 1.55001], [2.2799899999999997, 2.32001], [1.24999, 1.55001]]

    lb, ub, ulb, lub = Kodiak.minmax_diff(expression, jac_mat, bias, bounds)
    print("Kodiak computed enclosure: [{}, {}, {}, {}]".format(lb, ub, ulb, lub))

    assert lb < ulb
    assert ub > lub

    # default precious is -3
    assert ulb - lb < 1e-2
    assert ub - lub < 1e-2

#def test_atan():
#    'test atan functionality'

#    Kodiak.add_variable('x')
#    Kodiak.add_variable('y')

#    e = Kodiak.to_expression("atan(x / y)")

#    for y_bound in [[-1, -1e-6], [1e-6, 1]]:
#        bounds = [[1, 1], y_bound]

#        lb, ub, ulb, lub = Kodiak.minmax_diff(e, [0, 0], 0, bounds)
#        dif_lb = ulb - lb
#        dif_ub = ub - lub
#        print("Kodiak computed enclosure: [{}, {}] with tolerances {} and {}".format(lb, ub, dif_lb, dif_ub))

#    assert False

#test_atan()
