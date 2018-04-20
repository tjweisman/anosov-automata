import representation
import numpy as np
from kbmag_graphs import get_graph_from_fsa


freegp_param = 3.0

#this is the image of the INVERSES of the generators under the rep
representation = {
    "a":np.matrix([[1.0, freegp_param],
                    [0.0, 1.0]]),
    "b":np.matrix([[1.0, 0.0],
                    [freegp_param, 1.0]])
}

dim = 2

for generator, matrix in representation.items():
    representation[generator.upper()] = np.linalg.inv(matrix)

def random_geodesic_walk(graph, vertex, length,
                         matrix = np.identity(2)):
    if length > 0:
        next_vertex = np.random.choice(graph.neighbors_in(vertex))
        matrix = representation[graph.edge_label(next_vertex, vertex)] * matrix
        return random_geodesic_walk(graph, next_vertex, length - 1, matrix)
    else:
        u,s,vh = np.linalg.svd(matrix)
        return u
    

free2 = get_graph_from_fsa("kbmag_files/f2.geowa")
free2.delete_vertices(flatten(free2.level_sets()))
