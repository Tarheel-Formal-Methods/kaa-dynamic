import numpy as np
import sympy as sp

from operator import add
from functools import reduce

from kaa.parallelotope import Parallelotope
from kaa.lputil import minLinProg, maxLinProg
from kaa.settings import KaaSettings

from kaa.timer import Timer

OptProd = KaaSettings.OptProd

class Bundle:

    def __init__(self, model, T, L, offu, offl):

        assert np.size(L,0) == np.size(offu,0), "Directions matrix L and upper offsets must have matching dimensions: {} {}".format(np.size(L,0), np.size(offu,0))
        assert np.size(L,0) == np.size(offl,0), "Directions matrix L and lower offsets must have matching dimensions {} {}".format(np.size(L,0), np.size(offl,0))
        assert np.size(T,1) == np.size(L,1), "Template matrix T must have the same dimensions as Directions matrix L"

        self.T = T 
        self.L = L
        self.offu = offu
        self.offl = offl
        self.model = model
        self.vars = model.vars
        self.sys_dim = len(model.vars)
        self.num_direct = len(self.L)


    """
    Returns linear constraints representing the polytope defined by bundle.
    @returns linear constraints and their offsets.
    """
    def getIntersect(self):
        A = np.empty([2*self.num_direct, self.sys_dim])
        b = np.empty(2*self.num_direct) 

        for ind in range(self.num_direct):
            A[ind] = self.L[ind]
            A[ind + self.num_direct] = np.negative(self.L[ind])
            b[ind] = self.offu[ind]
            b[ind + self.num_direct] = self.offl[ind]

        return A, b

    """
    Returns the bundle with tightest offsets for each direction vector in self.L
    i.e each hyperplane defined by the direction vector is re-fitted to be tangent to the polytope.
    @returns canonized Bundle object
    """
    def canonize(self):
        A, b = self.getIntersect()

        for row_ind, row in enumerate(self.L):

            self.offu[row_ind] = (maxLinProg(row, A, b)).fun
            self.offl[row_ind] = (maxLinProg(np.negative(row), A, b)).fun

    """
    Returns the Parallelotope object defined by a row in the template matrix.
    @params temp_ind: index of row corresponding to desired parallelotope.i
    @returns Parallelotope object described by T[temp_ind]
    """
    def getParallelotope(self, temp_ind):

        A = np.empty([2*self.sys_dim,self.sys_dim])
        b = np.empty(2*self.sys_dim)

        'Fetch linear constraints defining the parallelotope.'
        for fac_ind, facet in enumerate(self.T[temp_ind].astype(int)):
            A[fac_ind] = self.L[facet]
            A[fac_ind + self.sys_dim] = np.negative(self.L[facet])
            b[fac_ind] = self.offu[facet]
            b[fac_ind + self.sys_dim] = self.offl[facet]

        return Parallelotope(A, b, self.vars)

    """
    Add a matrix of templates to the end of templates matrix.
    @params temp_row_mat: Matrix of new template entries.
    """
    def add_temp(self, temp_row_mat):
        self.T = np.append(self.T, temp_row_mat, axis=0)

    """
    Remove specified template entries from templates matrix.
    @params temp_idx: list of indices specifying row indices in template matrix
    """
    def remove_temp(self, temp_idx):
        self.T = np.delete(self.T, temp_idx, axis=0)

    """
    Add matrix of direction to end of directions matrix.
    Appends new row elemets into the offset matrices as well.
    @params dir_row_mat: Matrix of new direction entries
    """
    def add_dir(self, dir_row_mat):
        A, b = self.getIntersect()

        self.L = np.append(self.L, dir_row_mat, axis=0)
        new_uoffsets = [[ (maxLinProg(row, A, b)).fun for row in dir_row_mat ]]
        new_loffsets = [[ (maxLinProg(np.negative(row), A, b)).fun for row in dir_row_mat ]]
        
        self.offu = np.append(self.offu, new_uoffsets)
        self.offl = np.append(self.offl, new_loffsets)
        self.num_direct = len(self.L)
        
    """
    Remove specified direction entries from directions matrix.
    @params temp_idx: list of indices specifying row indices in directions matrix
    """
    def remove_dir(self, dir_idx):
        self.L = np.delete(self.L, dir_idx, axis=0)
        self.num_direct = len(self.L)

        self.offu = np.delete(self.offu, dir_idx, axis=0)
        self.offl = np.delete(self.offl, dir_idx, axis=0)

    """
    Copies bundle information into a new Bundle object
    @returns new Bundle object copy.
    """
    def copy(self):
        new_T = np.copy(self.T)
        new_L = np.copy(self.L)
        new_offu = np.copy(self.offu)
        new_offl = np.copy(self.offl)
        
        return Bundle(self.model, new_T, new_L, new_offu, new_offl)

"""
Object responsible for transforming Bundle objects according to input model's dynamics.
"""
class BundleTransformer:

    def __init__(self, model):
        self.f = model.f

    """
    Transforms the bundle according to the dynamics governing the system. (dictated by self.f)

    @params bund: Bundle object to be transformed under dynamics.
    @returns canonized transformed bundle.
    """
    def transform(self, bund):
        new_offu = np.full(bund.num_direct, np.inf)
        new_offl = np.full(bund.num_direct, np.inf)

        for row_ind, row in enumerate(bund.T):
            
            'Find the generator of the parallelotope.'
            p = bund.getParallelotope(row_ind)
            genFun = p.getGeneratorRep()

            'Create subsitutions tuples.'
            var_sub = []
            for var_ind, var in enumerate(bund.vars):
                var_sub.append((var, genFun[var_ind]))

            for column in row.astype(int):
                curr_L = bund.L[column]

                'Perform functional composition with exact transformation from unitbox to parallelotope.'
                Timer.start('Functional Composition')
                fog = [ f.subs(var_sub) for f in self.f ]
                Timer.stop('Functional Composition')

                bound_polyu = [ curr_L[func_ind] * func for func_ind, func in enumerate(fog) ]
                bound_polyu = reduce(add, bound_polyu)

                'Calculate min/max Bernstein coefficients.'
                Timer.start('Kodiak Computation')
                ub, lb = OptProd(bound_polyu, bund).getBounds()
                Timer.stop('Kodiak Computation')

                new_offu[column] = min(ub, new_offu[column])
                new_offl[column] = min(-1 * lb, new_offl[column])

        bund.offu = new_offu
        bund.offl = new_offl
        bund.canonize()

        return bund
