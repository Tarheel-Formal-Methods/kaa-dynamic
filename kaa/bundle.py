import numpy as np
import sympy as sp
from enum import Enum
import warnings

from kaa.parallelotope import Parallelotope
from kaa.templates import TempStrategy
from kaa.linearsystem import LinearSystem
from kaa.lputil import minLinProg, maxLinProg
from kaa.settings import KaaSettings
from kaa.timer import Timer

OptProd = KaaSettings.OptProd

warnings.filterwarnings('ignore')

class Bundle:

    def __init__(self, model, T, L, offu, offl):

        assert np.size(L,0) == np.size(offu,0), "Directions matrix L and upper offsets must have matching dimensions: {} {}".format(np.size(L,0), np.size(offu,0))
        assert np.size(L,0) == np.size(offl,0), "Directions matrix L and lower offsets must have matching dimensions {} {}".format(np.size(L,0), np.size(offl,0))
        assert np.size(T,1) == np.size(L,1), "Template matrix T must have the same dimensions as Directions matrix L"

        'Label the initial directions and templates with Default moniker.'
        self.labeled_L = [ (dir_row, "Default" + str(row_idx)) for row_idx, dir_row in enumerate(L) ]
        self.labeled_T = np.asarray([(row, "DefaultTemp", 1) for row in self.__convert_to_labeled_T(T)])

        self.offu = offu
        self.offl = offl

        self.model = model
        self.vars = model.vars
        self.dim = model.dim

        self.num_dir = len(L)
        self.num_temp = len(T)
        self.num_strat = 1

        self.strat_temp_id = {}

    @property
    def T(self):

        curr_T = self.__get_row(self.labeled_T)
        gen_T = np.empty((self.num_temp, self.dim))
        #print(f"L: {self.labeled_L}")

        for row_idx, row_labels in enumerate(curr_T):
            for col_idx, row_label in enumerate(row_labels):
                gen_T[row_idx][col_idx] = self.__get_dir_row_from_label(row_label)

        return gen_T

    @property
    def L(self):
        return np.asarray(self.__get_row(self.labeled_L))

    "Returns list of Parallelotope objects defining this bundle. WARNING: superfluous calls to getParallelotope for now"
    @property
    def ptopes(self):
        return [self.getParallelotope(i) for i in range(self.num_temp)]
    
    """
    Returns linear constraints representing the polytope defined by bundle.
    @returns linear constraints and their offsets.
    """
    def getIntersect(self):

        L = self.L
        A = np.empty([2*self.num_dir, self.dim])
        b = np.empty(2*self.num_dir)

        for ind in range(self.num_dir):
            A[ind] = L[ind]
            A[ind + self.num_dir] = np.negative(L[ind])
            b[ind] = self.offu[ind]
            b[ind + self.num_dir] = self.offl[ind]

        return LinearSystem(self.model, A, b)


    """
    Returns the bundle with tightest offsets for each direction vector in self.L
    i.e each hyperplane defined by the direction vector is re-fitted to be tangent to the polytope.
    @returns canonized Bundle object
    """
    def canonize(self):
        bund_sys = self.getIntersect()
        L = self.L

        for row_ind, row in enumerate(L):
            self.offu[row_ind] = bund_sys.max_opt(row).fun
            self.offl[row_ind] = bund_sys.max_opt(np.negative(row)).fun

    """
    Returns list of Parallelotopes by the strategy they are associated with.
    """
    def get_ptopes_by_strat(self, strat):
        assert isinstance(strat, TempStrategy), "input must be TempStrategy"

        if str(strat) not in self.strat_temp_id:
            return [None]

        temp_id = self.strat_temp_id[str(strat)]

        asso_temps = []
        for temp_idx, (_,_, tid) in enumerate(self.labeled_T):
            if tid == temp_id:
                asso_temps.append(temp_idx)

        return [self.getParallelotope(temp_idx) for temp_idx in asso_temps]


    """
    Returns the Parallelotope object defined by a row in the template matrix.
    @params temp_ind: index of row corresponding to desired parallelotope.i
    @returns Parallelotope object described by T[temp_ind]
    """
    def getParallelotope(self, temp_ind):
        L = self.L
        T = self.T
        A = np.empty([2*self.dim,self.dim])
        b = np.empty(2*self.dim)

        'Fetch linear constraints defining the parallelotope.'
        for fac_ind, facet in enumerate(T[temp_ind].astype(int)):
            A[fac_ind] = L[facet]
            A[fac_ind + self.dim] = np.negative(L[facet])
            b[fac_ind] = self.offu[facet]
            b[fac_ind + self.dim] = self.offl[facet]

        return Parallelotope(self.model, A, b)

    """
    Add a matrix of templates to the end of templates matrix.
    @params temp_row_mat: Matrix of new template entries.
    """
    def add_temp(self, asso_strat, row_labels, temp_label):
        assert len(row_labels) == self.dim, "Number of directions to use in template must match the dimension of the system."

        new_temp_ent = (self.__get_global_labels(asso_strat, row_labels), self.__get_global_labels(asso_strat, temp_label), self.__get_temp_id(asso_strat))
        self.labeled_T = np.append(self.labeled_T, [new_temp_ent], axis=0)
        self.num_temp = len(self.labeled_T)

    """
    Remove specified template entries from templates matrix.
    @params temp_idx: list of indices specifying row indices in template matrix
    """
    def remove_temp(self, asso_strat, temp_label):
        label_indices = self.__get_label_indices(self.labeled_T, self.__get_global_labels(asso_strat, temp_label))
        self.labeled_T = np.delete(self.labeled_T, label_indices, axis=0)
        #print(f"labeled_T: {self.labeled_T}")
        self.num_temp = len(self.labeled_T)

    """
    Add matrix of direction to end of directions matrix.
    Appends new row elemets into the offset matrices as well.
    @params dir_row_mat: Matrix of new direction entries

    """
    def add_dirs(self, asso_strat, dir_row_mat, dir_labels):
        assert len(dir_row_mat) == len(dir_labels), "Number of input direction rows must be one-to-one with the labels"
        bund_sys = self.getIntersect()
        prev_len = self.num_dir

        'Update new templates to envelope current polytope'
        labeled_L_ents = zip(dir_row_mat, self.__get_global_labels(asso_strat, dir_labels))
        self.labeled_L = np.append(self.labeled_L, list(labeled_L_ents), axis=0)

        new_uoffsets = [[ bund_sys.max_opt(row).fun for row in dir_row_mat ]]
        new_loffsets = [[ bund_sys.max_opt(np.negative(row)).fun for row in dir_row_mat ]]
        
        self.offu = np.append(self.offu, new_uoffsets)
        self.offl = np.append(self.offl, new_loffsets)
        self.num_dir = len(self.labeled_L)
        
    """
    Remove specified direction entries from directions matrix from their labels.
    @params temp_idx: list of indices specifying row indices in directions matrix
    """
    def remove_dirs(self, asso_strat, labels):
    
        label_indices = self.__get_label_indices(self.labeled_L, self.__get_global_labels(asso_strat, labels))

        #print(f"RemoveDir: Labels: {global_labels}")
        #print(f"labeled_L: {self.labeled_L}")
        #print(f"Mask: {label_indices}")
        
        self.labeled_L = np.delete(self.labeled_L, label_indices, axis=0)
        self.num_dir = len(self.labeled_L)

        self.offu = np.delete(self.offu, label_indices, axis=0)
        self.offl = np.delete(self.offl, label_indices, axis=0)


    def __get_temp_id(self, asso_strat):

        if str(asso_strat) not in self.strat_temp_id:
            self.num_strat += 1
            self.strat_temp_id[str(asso_strat)] = self.num_strat

        return self.strat_temp_id[str(asso_strat)]

    """
    Converts relative labels given by strategies to global labels understood by the Bundle object.
    All labels will generally be in the form of strategy_name + relative_label
    @params: asso_strat: strategy which claims the set of labels or label
             labels: set of labels or label which must be converted into the global format understood by Bundle
    @returns list of converted labels or label.
    """
    def __get_global_labels(self, asso_strat, labels):
        return [str(asso_strat) + label for label in labels] if isinstance(labels, list) else str(asso_strat) + labels

    """
    Returns the indices which a list of indices or an index points towrads in the tuple matrix.
    @params tup mat: input tuple matrix
            labels: list of labels or a single label relevant to search
    @returns list of indices which a label in the input points towards
    """
    def __get_label_indices(self, tup_mat, labels):
        label_set = set(labels) if isinstance(labels, list) else set([labels])
        #print(self.__get_label(tup_mat))
        return [idx for idx, label in enumerate(self.__get_label(tup_mat)) if label in label_set]

    """
    Returns the row data in order of the labeled matrix rows. The input will be a
    labeled matrix i.e list of tuples where the first coordinate contains the data
    and the second value contains the label of the data.
    @params tup_mat: input tuple matrix i.e self.labeled_L or self.labeled_T
    @returns list of rows (data) in order which they appear in the input tuple matrix.
    """
    def __get_row(self, tup_mat):
        return list(map(lambda row_label_tup: row_label_tup[0], tup_mat))

    """
    Returns the labels in order of the labeled matrix rows. The input will be a
    labeled matrix i.e list of tuples where the first coordinate contains the data
    and the second value contains the label of the data.
    @params tup_mat: input tuple matrix i.e self.labeled_L or self.labeled_T
    @returns list of labels in order which they appear in the input tuple matrix.
    """
    def __get_label(self, tup_mat):
        return list(map(lambda row_label_tup: row_label_tup[1], tup_mat))

    """
    Searches for the row which the input label points to in the directions matrix.
    @params dir_label: label to search for
    @returns index which label points to.
    """
    def __get_dir_row_from_label(self, dir_label):
        label_indices = self.__get_label_indices(self.labeled_L, dir_label)
        assert len(label_indices) == 1, f"Every index should have a unique label. \n Label: {dir_label} \n Offending list: {label_indices}  \n Bundle label dump: \n {self.labeled_L}"
        return label_indices[0]

    """
    Converts template matrix from numpy form into a labeled form
    containing only the list of labels pointing to the relevant rows in
    the directions matrix
    @params T: numpy template matrix with integer entries
    @returns list of lists containing labels pointing to direction matrix rows.
    """
    def __convert_to_labeled_T(self, T):
        labeled_T = []
        for temp_row in T:
            labeled_T.append([self.labeled_L[col][1] for col in temp_row.astype(int)])

        return np.asarray(labeled_T)

    """
    Copies bundle information into a new Bundle object
    @returns new Bundle object copy.

    def copy(self):
        new_T = np.copy(self.T)
        new_L = np.copy(self.L)
        new_offu = np.copy(self.offu)
        new_offl = np.copy(self.offl)
        
        return Bundle(self.model, new_T, new_L, new_offu, new_offl)
    """
    
"""
Wrapper over bundle transformation modes.
"""
class BundleMode(Enum):
    AFO = False
    OFO = True

"""
Object responsible for transforming Bundle objects according to input model's dynamics.
"""
class BundleTransformer:

    """
    Constuctor for BundleTransformer
    @params mode: mode for performing bundle transformations.
    A value of BundleMode.OFO (1) indicates using the One-for-One transformation method.
    Otherwise, a value of BundleMode.AFO (0) indicates using the All-for-One transformation method.
    """
    def __init__(self, model, mode):
        self.f = model.f
        self.vars = model.vars
        self.ofo_mode = mode

    """
    Transforms the bundle according to the dynamics governing the system. (dictated by self.f)

    @params bund: Bundle object to be transformed under dynamics.
    @returns canonized transformed bundle.
    """
    def transform(self, bund):
        L = bund.L
        T = bund.T
        new_offu = np.full(bund.num_dir, np.inf)
        new_offl = np.full(bund.num_dir, np.inf)

        for row_ind, row in enumerate(T):
            'Find the generator of the parallelotope.'
            ptope = bund.getParallelotope(row_ind)
            #print(f"Ptope {row_ind}\n")

            'Toggle iterators between OFO/AFO'
            direct_iter = row.astype(int) if self.ofo_mode.value else range(bund.num_dir)
            #print(f"Num of directions: {bund.num_dir}")

            for column in direct_iter:
                curr_L = L[column]
                ub, lb = self.__find_bounds(curr_L, ptope, bund)

                new_offu[column] = min(ub, new_offu[column])
                new_offl[column] = min(lb, new_offl[column])

        bund.offu = new_offu
        bund.offl = new_offl

        #print("Upper Offsets: {}\n".format(bund.offu))
        #print("Lower Offsets: {}\n".format(bund.offl))

        bund.canonize()
        return bund

    """
    Find bounds for max c^Tf(x) over paralleltope
    @params: ptope: Parallelotope object to optimize over.
             dir_vec: direction vector
    @returns: upper bound, lower bound
    """
    def __find_bounds(self, dir_vec, ptope, bund):

        'Find the generator of the parallelotope.'
        genFun = ptope.getGeneratorRep()

        'Create subsitutions tuples.'
        var_sub = [(var, genFun[var_ind]) for var_ind, var in enumerate(self.vars)]

        #print(f"Variable Sub for {dir_vec}: {var_sub}")
        
        Timer.start('Functional Composition')
        fog = [func.subs(var_sub, simultaneous=True) for func in self.f]
        Timer.stop('Functional Composition')

        'Perform functional composition with exact transformation from unitbox to parallelotope.'
        bound_polyu = 0
        for coeff_idx, coeff in enumerate(dir_vec):
            bound_polyu += coeff * fog[coeff_idx]

        'Calculate min/max Bernstein coefficients.'
        Timer.start('Bound Computation')
        ub, lb = OptProd(bound_polyu, bund).getBounds()
        Timer.stop('Bound Computation')

        return ub, -1 * lb
