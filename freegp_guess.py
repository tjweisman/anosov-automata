import numpy as np
import matplotlib.pyplot as plt

import pyximport; pyximport.install()
import multicone_guess

from kbmag_utils import load_fsa_dict
from representation import Representation
from fsa import FSA

def freegp_rep():
    free2 = FSA(load_fsa_dict("kbmag_files/f2.geowa"))
    free2.make_recurrent()

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
    plt.show()
