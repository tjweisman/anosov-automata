import numpy as np

from convexhull import ConvexHull
from projective_geom import *


class Multicone:
    def __init__(self, pts):
        self._pts = pts
        self._affine_pts, self._chart = affine_coords(pts)

        self.recompute_hull()

    def contains(self, pts):
        vertices = affine_coords(pts, self._chart)
        limiting_matrix = np.column_stack([vertices,
                                           np.ones(len(vertices))])
        
        return (np.matmul(self._hull.equations, limiting_matrix.T) < 0).all()

    def expand_to_include(self, multicone):
        pts = np.concatenate([self._pts, multicone.points], axis=0)
        return Multicone(pts)

    def expand_by_epsilon(self, epsilon):
        #TODO: make this expansion by an actual metric on RP^n
        center = np.mean(self._affine_pts, axis=0)
        affine_pts = (self._affine_pts - center) * (1 + epsilon) + center
        pts = lift(affine_pts, self._chart)
        return Multicone(pts)


    def recompute_hull(self):
        hull = ConvexHull(self._affine_pts)
        self._pts = self._pts[hull.vertices]
        self._affine_pts = self._affine_pts[hull.vertices]
        self._hull = ConvexHull(self._affine_pts)

    def transform(self, matrix):
        return Multicone(self._pts * matrix.T)

    @property
    def points(self):
        return self._pts
