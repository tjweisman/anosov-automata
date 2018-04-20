import numpy as np
import matplotlib.pyplot as plt

import pyximport; pyximport.install()
import multicone_guess

from kbmag_graphs import get_graph_from_fsa
from representation import Representation

def freegp_rep():
    free2 = get_graph_from_fsa("kbmag_files/f2.geowa")
    free2.delete_vertices(flatten(free2.level_sets()))

    freegp_param = 3.0
    rho = Representation(["a", "b"])
    rho["a"] = np.matrix([[1.0, freegp_param],
                          [0.0, 1.0]])
    rho["b"] = np.matrix([[1.0, 0.0],
                          [freegp_param, 1.0]])

    return (rho, free2)

def plot_clusters(clusters, filename):
    color = iter(plt.cm.rainbow(np.linspace(0, 1, len(clusters))))
    for v, cluster in clusters.items():
        plt.scatter((cluster.T)[0], (cluster.T)[1], c=next(color))
    plt.savefig(filename)
