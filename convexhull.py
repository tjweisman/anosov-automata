import scipy.spatial
import numpy as np

class ConvexHull(scipy.spatial.ConvexHull):
    def __init__(self, pts):
        self.d1 = False
        if len(pts.shape) == 1 or pts.shape[1] == 1:
            self.d1 = True
            self._points = pts
            self.calc_1d_limits()
        else:
            super(ConvexHull, self).__init__(pts)

    def calc_1d_limits(self):
        if self.d1:
            pts = self._points
            self._vertices = [np.argmin(pts), np.argmax(pts)]
            self.volume = pts[self._vertices[1]][0] - pts[self._vertices[0]][0]
            self.equations = np.array([
                [1.0, -1 * pts[self._vertices[1]][0]],
                [-1.0, pts[self.vertices[0]][0]]
            ])
            
            
    def contains(self, hull):
        limiting_matrix = np.column_stack([hull.points[hull.vertices],
                                           np.ones(len(hull.vertices))])
        
        return (np.matmul(self.equations, limiting_matrix.T) < 0).all()
