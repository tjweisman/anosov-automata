import numpy as np
import random
from scipy.spatial import ConvexHull

from multicone import Multicone

def random_geodesic_walk(representation, graph, vertex, length):
    """singular directions for the endpoint of a long geodesic.

    Walk BACKWARDS randomly along the geodesic automaton, computing
    the image of the word spelled out under the representation. Then
    return singular directions of the resulting matrix.

    """
    if matrix is None:
        matrix = np.identity(representation.dim)
        
    for i in range(length):
        edges = list(graph.edges_in(vertex))
        next_vertex, vertex, label = random.choice(edges)
        step = representation.invert_gen(label)
        matrix = matrix * representation[step]
        vertex = next_vertex
        
    u,s,vh = np.linalg.svd(matrix)
    return u

def point_cluster(representation, graph, vertex, size, length=100):
    """get a cluster of singular values for long geodesics.

    Do a bunch of random walks ending at a particular vertex, and find
    a cluster of largest s.v.'s for these matrices

    """
    pts = np.zeros((size, representation.dim))

    for i in range(size):
        spaces = random_geodesic_walk(representation, graph, vertex, length)
        pts[i] = (spaces.T)[0]

    return pts

def build_clusters(representation, graph, size, length=100):
    """get clusters of points (hopefully approximating multicones)

    for each vertex of the automaton, get a cluster of largest
    singular directions coming from long (backwards) walks on the
    automaton ending at that vertex.

    """
    clusters = {}
    for v in graph.vertices():
        clusters[v] = point_cluster(representation, graph, v, size, length)
    return clusters

def guess_multicones(rep, graph):
    """get an approximation for a multicone over each vertex of the
    automaton.

    

    """
    clusters = build_clusters(rep, graph, 500)
    
    multicone_guesses = {v:Multicone(cluster)
                         for v, cluster in clusters.iteritems()}
    
    return multicone_guesses

def expand_multicones(rep, graph, cones, epsilon = 0.01, iterations = 3):
    expanded_cones = {v:cone for v,cone in cones.iteritems()}
    for i in range(iterations):
        for v in graph.vertices():
            for vt, w, label in graph.edges_out(v):
                transformed = (
                    expanded_cones[v].transform(rep[rep.invert_gen(label)])
                )
                if not expanded_cones[w].contains(transformed.points):
                    expanded_cones[w] = (
                        expanded_cones[w].expand_to_include(transformed)
                    )
                    expanded_cones[w] = (
                        expanded_cones[w].expand_by_epsilon(epsilon)
                    )
                    
    return expanded_cones


def check_inclusions(rep, graph, multicones):
    for v in graph.vertices():
        for vt, w, label in graph.edges_out(v):
            v_pts = multicones[v].points
            mat = rep[rep.invert_gen(label)]
            transformed_cone = v_pts * mat.T
            
            if not multicones[w].contains(transformed_cone):
                print("edge {0} --{2}--> {1} does not "
                      "satisfy inclusion".format(
                          vt, w, label))
