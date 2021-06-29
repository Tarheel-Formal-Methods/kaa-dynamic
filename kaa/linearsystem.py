import random as rand
import numpy as np
import multiprocessing as mp
from scipy.spatial import HalfspaceIntersection, ConvexHull
from scipy.spatial.qhull import QhullError
from operator import mul, add
from functools import reduce
from collections import namedtuple

from kaa.lputil import minLinProg, maxLinProg, LPUtil
from settings import KaaSettings
from kaa.trajectory import Traj, TrajCollection
from settings import KaaSettings
from kaa.log import Output

VolDataTuple = namedtuple('VolDataTuple', ['ConvHullVol', 'EnvelopBoxVol'])


class ChebyCenter:

    def __init__(self, center, radius):
        self.center = center
        self.radius = radius


class LinearSystem:

    def __init__(self, model, A, b, constr_mat=None):
        self.A = A
        self.b = b
        self.model = model
        self.vars = model.vars
        self.dim = model.dim
        self.constr_mat = constr_mat  # Pointer to total constraint mat for LPUtil purposes.
        self.randgen = rand.Random(KaaSettings.RandSeed)

    """
    Computes and returns the Chebyshev center of parallelotope.
    @returns self.dim point marking the Chebyshev center.
    """
    @property
    def chebyshev_center(self):
        'Initialize objective function for Chebyshev intersection LP routine.'
        c = [0 for _ in range(self.dim + 1)]
        c[-1] = 1

        row_norm = np.reshape(np.linalg.norm(self.A, axis=1), (self.A.shape[0], 1))
        center_A = np.hstack((self.A, row_norm))

        center_pt = maxLinProg(self.model, c, center_A, self.b).x
        return ChebyCenter(center_pt[:-1], center_pt[-1])

    """
    Volume estimation of system by sampling points and taking ratio.
    @params samples: number of samples used to estimate volume
    @returns estimated volume of linear system stored in VolDataTuple
    """
    @property
    def volume(self):
        envelop_box_vol = self.calc_vol_envelop_box()
        conv_hull_vol = self.calc_vol_conv_hull() if self.dim < 4 else None

        return VolDataTuple(conv_hull_vol, envelop_box_vol)

    """
    Find vertices of this linear system.
    """
    @property
    def vertices(self):
        phase_intersect = np.hstack((self.A, - np.asarray([self.b]).T))

        'Run scipy.spatial.HalfspaceIntersection.'
        try:
            center_pt = np.asarray(self.chebyshev_center.center)
            hs = HalfspaceIntersection(phase_intersect, center_pt)
        except QhullError:
            feasible_pt = self.feasible_point()
            hs = HalfspaceIntersection(phase_intersect, feasible_pt)

        vertices = np.asarray(hs.intersections)
        return vertices

    def feasible_point(self, use_interior_point=True, num_samples=5):
        if use_interior_point:
            return minLinProg(self, np.zeros(self.dim), self.A, self.b, method='Interior').x
        else:
            sample_mat = np.random.randn(num_samples, self.dim)

            perturbed_feasible_pts = np.empty((num_samples, self.dim))
            for idx, perturbed_obj in enumerate(sample_mat):
                perturbed_feasible_pts[idx] = self.min_obj(perturbed_obj).x

            return np.mean(perturbed_feasible_pts, axis=0)

    """
    Calculate the volume of the smallest enveloping box of linear system.
    """
    def calc_vol_envelop_box(self):
        envelop_box = self.calc_envelop_box()
        return calc_box_volume(envelop_box)

    """
    Calculate the volume of the convex hull of linear system.
    """
    def calc_vol_conv_hull(self):
        return ConvexHull(self.vertices).volume

    """
    Maxmize optimization function y over Ax \leq b
    @params y: linear function to optimize over
    @returns LinProgResult
    """
    def max_obj(self, y):
        assert len(y) == self.dim, "Linear optimization function must be of same dimension as system."
        return maxLinProg(self.model, y, self.A, self.b, constr_mat=self.constr_mat)

    """
    Minimize optimization function y over Ax \leq b
    @params y: linear function to optimize over
    @returns LinProgResult
    """
    def min_obj(self, y):
        assert len(y) == self.dim, "Linear optimization function must be of same dimension as system."
        return minLinProg(self.model, y, self.A, self.b, constr_mat=self.constr_mat)

    """
    Checks if point is indeed contained in Ax \leq b
    @params point: point to test
    @returns boolean value indictating membership.
    """
    def check_membership(self, point):
        assert len(point) == self.dim, "Point must be of the same dimension as system."
        for row_idx, row in enumerate(self.A):
            if np.dot(row, point) > self.b[row_idx]:
                return False
        return True

    """
    Calculate the enveloping box over the linear system
    @params model: input model
    @returns list of intervals representing edges of box.
    """
    def calc_envelop_box(self):
        box_interval = np.empty((self.dim,2))

        for i in range(self.dim):
            y = np.zeros(self.dim)
            y[i] = 1

            maxCood = self.max_obj(y).fun
            minCood = self.min_obj(y).fun
            box_interval[i] = [minCood, maxCood]

        return box_interval

    """
    Sample a random point contained within a box.
    @params box_intervals: list of lists defining box
    @returns random point sampled within box
    """
    def __sample_box_points(self, box_intervals, num_points):
        assert len(box_intervals) == self.dim, "Number of intervals defining box must match dimension of system."
        return [[self.randgen.uniform(start, end) for start, end in box_intervals] for _ in range(num_points)]

    """
    Create a trajectory at certain point in system and propagate it forward according
    to system dynamics. Add to shared Queue after initialization is finished.
    @params model: system model
            point: initial point of trajectory
            time_step: number of steps to propagate forward
            output_queue: Manager.Queue object shared between processes
    """
    def create_traj(self, point, steps, output_queue=None):
        if output_queue:
            output_queue.put(Traj(self.model, point, steps=steps))
        else:
            return Traj(self.model, point, steps=steps)

    """
    Auxiliary method to convert Manager.Queue object into list.
    @params mp_queue: Manager.Queue object to convert.
    @returns: List object with queue contents.
    """
    def queue_to_list(self, mp_queue):
        output_list = []
        while(mp_queue.qsize() != 0):
            output_list.append(mp_queue.get())

        return output_list

    """
    Generate random trajectories from polytope defined by parallelotope bundle.
    @params model: Model
            num_traj: number of trajectories to generate.
            time_steps: number of time steps to generate trajs.
    @returns list of Traj objects representing num random trajectories.
    """
    def generate_ran_trajs(self, num_trajs, steps):
        output_trajs = TrajCollection(self.model)
        initial_points = self.sample_ran_pts_envelop_box(num_trajs)

        if KaaSettings.Parallelize:
            'Parallelize point propagation'

            output_queue = mp.Manager().Queue()
            input_params = [(point, steps, output_queue) for point in initial_points]

            p = mp.Pool()
            p.starmap(self.create_traj, input_params)
            p.close()
            p.join()

            output_trajs.add(*self.queue_to_list(output_queue))
        else:
            for point in initial_points:
                output_trajs.add(self.create_traj(point, steps))

        return output_trajs

    """
    Worker for trajectory generation using support points.
    """
    def generate_supp_worker(self, dir_vec, steps, output_queue=None):
        supp_point = self.max_obj(dir_vec).x
        neg_supp_point = self.min_obj(dir_vec).x

        return [self.create_traj(supp_point, steps, output_queue),
                self.create_traj(neg_supp_point, steps, output_queue)]

    """
    Generates trajectories based on support points of provided directions.
    @params dir_vecs: Matrix of directions stored as rows.
            time_steps: number of steps to
    """
    def generate_supp_trajs(self, dir_vecs, steps):
        output_trajs = TrajCollection(self.model)
        if KaaSettings.Parallelize:
            output_queue = mp.Manager().Queue()

            input_params = []
            for dir_vec in dir_vecs:
                input_params += [(dir_vec, steps, output_queue), (np.negative(dir_vec), steps, output_queue)]

            p = mp.Pool()
            p.starmap(self.generate_supp_worker, input_params)
            p.close()
            p.join()

            output_trajs.add(*self.queue_to_list(output_queue))
        else:
            'Exploiting Warm-start LP'
            with LPUtil(self.model, None, self.A, self.b, None, 'Simplex') as lp_inst:
                lp_inst.populate_consts()
                for dir_vec in dir_vecs:
                    lp_inst.c = dir_vec
                    lp_inst.populate_obj_vec()

                    supp_point = lp_inst.solve("Max").x
                    neg_supp_point = lp_inst.solve("Min").x

                    output_trajs.add(self.create_traj(supp_point, steps),
                                     self.create_traj(neg_supp_point, steps))

        return output_trajs

    def sample_ran_pts_envelop_box(self, num_trajs):
        box_intervals = self.__calc_envelop_box()
        return self.__sample_box_points(box_intervals, num_trajs)

    """
    Generates random points contained within the tightest enveloping parallelotope of the Chevyshev sphere.
    @params bund: Bundle object
            num_trajs: number of trajs to generate
            shrinkfactor: factor to shrink the radius of the sphere. This allows a box with smaller dimensions
    @returns list of generated random points.
    """
    def sample_ran_pts_centered_box(self, num_trajs, shrinkfactor=1):
        chebycenter = self.chebyshev_center

        center = chebycenter.center
        radius = chebycenter.radius
        box_intervals = [[c - (radius*shrinkfactor), c + (radius*shrinkfactor)] for c in center]

        return self.__sample_box_points(box_intervals, num_trajs)




"""
Calculate the volume of the box described by a list of tuples.
@params box_intervals: List of tuples where the first component is the lower end
                       and the second component is the upper end.
@returns volume of box
"""
def calc_box_volume(box_intervals):
    box_dim = [abs(end - start) for start,end in box_intervals]
    return reduce(mul, box_dim)

def intersect(model, *lin_sys):
    assert all(isinstance(sys, LinearSystem) for sys in lin_sys), "Intersection operation only viable with another LinearSystem object."

    all_A = [sys.A for sys in lin_sys]
    all_b = [sys.b for sys in lin_sys]

    new_A = np.concatenate(all_A, axis=0)
    new_b = np.concatenate(all_b, axis=0)

    return LinearSystem(model, new_A, new_b)
