'''
Wrapper interface to GLPK
'''
import swiglpk as glpk
import numpy as np
from scipy.optimize import linprog
from itertools import product

from kaa.log import Output
from kaa.timer import Timer

class LPSolution:

    def __init__(self, x, fun):
        self.x = x
        self.fun = fun

minLinProg = lambda c, A, b: __swiglpk_linprog(c, A, b, "min")
maxLinProg = lambda c, A, b: __swiglpk_linprog(c, A, b, "max")

def __scipy_linprog(c, A, b, obj):
    if obj == "min":
        obj_c_vec = c
    elif obj == "max":
        obj_c_vec = np.negative(c)
    else:
        raise RuntimeException("Obj should either be min or max")

    print(obj_c_vec)
    lin_result = linprog(obj_c_vec, A_ub=A, b_ub=b, bounds=(None,None), method='revised simplex')
    return LPSolution(lin_result.x, lin_result.fun)


def __swiglpk_linprog(c, A, b, obj):
    #Output.bold_write(f"Started LP with c: {c}")
    Timer.start("LP Routines")
    glp_obj = glpk.GLP_MIN if obj == "min" else glpk.GLP_MAX

    lp = glpk.glp_create_prob()
    glpk.glp_set_obj_dir(lp, glp_obj)

    params = glpk.glp_smcp()
    glpk.glp_init_smcp(params)
    params.msg_lev = glpk.GLP_MSG_ERR

    #params.tm_lim = int(glpk.GLPK_TIMEOUT * 1000)
    params.out_dly = 2 * 1000 # start printing to terminal delay
    params.meth = glpk.GLP_DUAL

    num_rows = A.shape[0]
    num_cols = A.shape[1]
    mat_size = num_rows * num_cols

    glpk.glp_add_rows(lp, num_rows)

    for row_ind in range(num_rows):
        glpk.glp_set_row_bnds(lp, row_ind+1, glpk.GLP_UP, 0.0, float(b[row_ind]))

    glpk.glp_add_cols(lp, num_cols)

    for col_ind in range(num_cols):
        glpk.glp_set_col_bnds(lp, col_ind+1, glpk.GLP_FR, 0.0, 0.0)
        glpk.glp_set_obj_coef(lp, col_ind+1, c[col_ind])

    'Swig arrays are used for feeding constraints in GLPK'

    ia, ja, ar = [],[],[]
    for i,j in product(range(num_rows),range(num_cols)):
        ia.append(i+1)
        ja.append(j+1)
        ar.append(float(A[i][j]))

    ia = glpk.as_intArray(ia)
    ja = glpk.as_intArray(ja)
    ar = glpk.as_doubleArray(ar)

    glpk.glp_load_matrix(lp, mat_size, ia, ja, ar)
    glpk.glp_simplex(lp, params)

    fun = glpk.glp_get_obj_val(lp)
    x = [i for i in map(lambda x: glpk.glp_get_col_prim(lp, x+1), range(num_cols))]

    glpk.glp_delete_prob(lp)
    glpk.glp_free_env()

    Timer.stop("LP Routines")
    #Output.bold_write(f"Ended LP with c: {c}")
    return LPSolution(x, fun)
