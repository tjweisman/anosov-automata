import numpy as np
from scipy.spatial import ConvexHull

def random_geodesic_walk(representation, graph, vertex, length,
                         matrix = None):
    if matrix is None:
        matrix = np.identity(representation.dim)
        
    for i in range(length):
        next_vertex = np.random.choice(graph.neighbors_in(vertex))
        step = representation.invert_gen(graph.edge_label(next_vertex, vertex))
        matrix = matrix * representation[step]
        vertex = next_vertex
        
    u,s,vh = np.linalg.svd(matrix)
    return u

def point_cluster(representation, graph, vertex, size, length=100):
    pts = np.zeros((size, representation.dim))

    for i in range(size):
        spaces = random_geodesic_walk(representation, graph, vertex, length)
        pts[i] = (spaces.T)[0]

    return pts

def build_clusters(representation, graph, size, length=100):
    clusters = {}
    for v in graph.vertices():
        clusters[v] = point_cluster(representation, graph, v, size, length)
    return clusters

def augment_clusters(clusters, representation, graph, iterations = 1):
    for x in range(iterations):
        for v in graph.vertices():
            for nbr in graph.neighbors_out(v):
                transform = representation[graph.edge_label(v, nbr)]
                np.append(clusters[nbr], (transform * clusters[v].T).T)

def compute_convex_hull(cluster):
    pass
