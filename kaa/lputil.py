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

'LP routines for problems without warm start'
minLinProg = lambda model, c, A, b, constr_mat: OneTimeLinProg(model, c, A, b, constr_mat, "Min")
maxLinProg = lambda model, c, A, b, constr_mat: OneTimeLinProg(model, c, A, b, constr_mat, "Max")

def OneTimeLinProg(model, c, A, b, constr_mat, glb_obj):
    with LPUtil(model, c, A, b, constr_mat) as lp_inst:
        lp_inst.populate_consts()
        lp_inst.populate_obj_vec()
        lp_sol = lp_inst.solve(glb_obj)

    return lp_sol

class LPUtil:

    def __init__(self, model, c=None, A=None, b=None, constr_mat=None):
        self.model = model
        self.dim = model.dim
        self.c = c
        self.A = A
        self.b = b
        self.constr_mat = constr_mat
        self.num_cols = None
        self.lp_prob = glpk.glp_create_prob()

        self.params = glpk.glp_smcp()
        glpk.glp_init_smcp(self.params)
        self.params.msg_lev = glpk.GLP_MSG_ERR
        self.params.meth = glpk.GLP_DUAL

        self.__normalize_constrs()

    def populate_consts(self):
        num_rows = self.A.shape[0]
        self.num_cols =  num_cols = self.A.shape[1]
        mat_size = num_rows * num_cols

        glpk.glp_add_rows(self.lp_prob, num_rows)

        for row_ind in range(num_rows):
            glpk.glp_set_row_bnds(self.lp_prob , row_ind+1, glpk.GLP_UP, 0.0, float(self.b[row_ind]))

        glpk.glp_add_cols(self.lp_prob , num_cols)

        'Swig arrays are used for feeding constraints in GLPK'
        ia, ja, ar = [],[],[]
        for i,j in product(range(num_rows),range(num_cols)):
            ia.append(i+1)
            ja.append(j+1)
            ar.append(float(self.A[i][j]))

        ia = glpk.as_intArray(ia)
        ja = glpk.as_intArray(ja)
        ar = glpk.as_doubleArray(ar)

        glpk.glp_load_matrix(self.lp_prob, mat_size, ia, ja, ar)

    def populate_obj_vec(self):
        for col_ind in range(self.num_cols):
            glpk.glp_set_col_bnds(self.lp_prob, col_ind+1, glpk.GLP_FR, 0.0, 0.0)
            glpk.glp_set_obj_coef(self.lp_prob, col_ind+1, self.c[col_ind])

    def solve(self, glp_obj):
        glp_obj = glpk.GLP_MAX if glp_obj == "Max" else glpk.GLP_MIN

        Timer.start("LP Routines")
        glpk.glp_set_obj_dir(self.lp_prob, glp_obj)

        glpk.glp_simplex(self.lp_prob, self.params)

        fun = glpk.glp_get_obj_val(self.lp_prob)
        x = [i for i in map(lambda x: glpk.glp_get_col_prim(self.lp_prob, x+1), range(self.num_cols))]

        Timer.stop("LP Routines")
        return LPSolution(x, fun)

    def __normalize_constrs(self):
        if self.constr_mat is not None:
            constr_norms = np.linalg.norm(self.A, axis=1)
            self.constr_mat /= constr_norms[:,None]

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        glpk.glp_delete_prob(self.lp_prob)
        glpk.glp_free_env()

"""
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
"""
