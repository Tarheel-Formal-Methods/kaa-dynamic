import numpy as np
import multiprocessing as mp
import random

from kaa.linearsystem import LinearSystem
from kaa.lputil import minLinProg, maxLinProg
from kaa.timer import Timer
from kaa.settings import KaaSettings

"""
Object encapsulating routines calculating properties of parallelotopes.
"""
class Parallelotope(LinearSystem):

    def __init__(self, model, A, b):
        super().__init__(model, A, b)
        self.u_A = A[:self.dim]
        self.u_b = b[:self.dim]

    @property
    def generator_vecs(self):
        bv, gens = self.__get_generators()
        return gens

    def __get_generators(self):
        base_vertex = self.__computeBaseVertex()
        return base_vertex, self.__computeGenerators(base_vertex)

    """
    Return list of functions transforming the n-unit-box over the parallelotope.
    @returns list of transfomation from unitbox over the parallelotope.
    """
    def getGeneratorRep(self):
        Timer.start('Generator Procedure')
        base_vertex, gens = self.__get_generators()

        'Create list representing the linear transformation q + \sum_{j} a_j* g_j where g_j'
        expr_list = base_vertex
        for j in range(self.dim):
            for var_ind, var in enumerate(self.vars):
                expr_list[j] += gens[var_ind][j] * var
        Timer.stop('Generator Procedure')

        return expr_list

    """
    Calculate generators as substraction: vertices - base_vertex.
    We calculate the vertices by solving the following linear system for each vertex i:


    Ax = [b_1, ... , -b_{i+n}, ... , b_n]^T

    Note that this simply finds the vertex to calculate the generator vectors. The generators will be the vectors g_j = v_j - q
    where v_j is the jth vertex calculated in the manner above.

    The parallelotope will be exprssed as sum of the base vertex and the convex combination of the generators.

    p(a_1, ... ,a_n) =  q + \sum_{j} a_j * g_j

    where q is the base vertex and the g_j are the generators. a_j will be in the unitbox [0,1]

    @params base_vertex: base vertex q
    @returns generator vectors g_j
    """
    def __computeGenerators(self, base_vertex):

        """
        'Hacky way to toggle parallelism for experiments'
        if KaaSettings.use_parallel:
            p = mp.Pool(processes=4)
            vertices = p.starmap(self._gen_worker, [ (i, u_b, coeff_mat) for i in range(self.dim) ])
            p.close()
            p.join()
        else:
        """
        vertices = []
        for i in range(self.dim):
            vertices.append(self.__gen_worker(i, self.u_b, self.u_A))

        vertex_list = [ [vert - base for vert, base in zip(vertices[i], base_vertex)] for i in range(self.dim) ]
        #print("Vertex List For Paratope: {} \n".format(vertices))
        #print("Vector List For Paratope: {} \n".format(vertex_list))
        return vertex_list

    """
    Worker process for calculating vertices of higher-dimensional parallelotopes.
    Only called by Pool.starmap
    @params i - vertex index
            u_b, coef_mat - shared reference to upper offsets and directions matrix.
    @returns coordinates of vertex
    """
    def __gen_worker(self, i, u_b, coeff_mat):
        #print(coeff_mat, u_b)

        negated_bi = np.copy(self.u_b)
        negated_bi[i] = -self.b[i + self.dim]
        #print("(coeff_mat, negated_bi): {}".format((coeff_mat, negated_bi)))
        sol_set_i = np.linalg.solve(coeff_mat, negated_bi)

        return list(sol_set_i)

    """
    Calculate the base vertex of the parallelotope (variable q)
    We calculate the vertices by solving a linear system of the following form:

    Ax = u_b

    where u_b are the offsets for the upper facets of parallelotope (first half of self.b).

    @returns base-vertex in list
    """
    def __computeBaseVertex(self):
        sol_set = np.linalg.solve(self.u_A, self.u_b)
        return list(sol_set)

    """
    Convert numpy matrix into sympy matrix
    @params mat: numpy matrix
    @returns sympy matrix counterpart
    """
    def __convertMatFormat(self, mat):
        return sp.Matrix(mat.tolist())

    """
    Takes solution set returned by sympy and converts into list
    @params fin_set: FiniteSet
    @returns list of sympy solution set
    """
    def __convertSolSetToList(self, fin_set):
        assert fin_set is not sp.EmptySet
        return list(fin_set.args[0])
