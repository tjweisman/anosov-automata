import numpy as np
from scipy.special import binom
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge


import pyximport; pyximport.install()
import multicone_builder

from multicone import Multicone
from projective_geom import sphere_lift
from kbmag_utils import load_fsa_dict
from representation import Representation
from fsa import FSA

def so_to_psl(A):
    #the isomorphism SO(2,1) --> PSL(2), assuming the matrix A is a
    #3x3 matrix determining a linear map in a basis where the
    #symmetric bilinear form has matrix diag(-1, -1, 1).
    
    a = np.sqrt((A[1][1] + A[2][1] + A[1][2] + A[2][2]) / 2)
    b = (A[0][1] + A[0][2]) / (2 * a)
    c = (A[1][0] + A[2][0]) / (2 * a)
    d = (A[2][0] - A[1][0]) / (2 * b)

    return np.matrix([[a, c],
                      [b, d]])

def psl_irrep(A, dim):
    #the irreducible representation from SL(2) to SL(dim) (via action
    #on homogeneous polynomials)
    
    a = A[0,0]
    b = A[0,1]
    c = A[1,0]
    d = A[1,1]
    
    im = np.matrix(np.zeros((dim, dim)))
    n = dim - 1
    for k in range(dim):
        for j in range(dim):
            for i in range(max(0, j - n + k), min(j+1, k+1)):
                im[j,k] += (binom(k,i) * binom(n - k, j - i)
                          * a**i * c**(k - i) * b**(j - i)
                          * d**(n - k - j + i))
    return im
    


def freegp_rep(trans_len = 50):
    #return FSA for a free group on 2 generators, together with a rep
    #consisting of a pair of matrices in SL(2) with orthogonal axes
    #and large translation length
    free2 = FSA(load_fsa_dict("kbmag_files/f2.geowa"))
    free2.make_recurrent()

    rtrans_len = 1.0 / trans_len
    
    rep = Representation(["a", "b"])
    rep["a"] = np.matrix([[trans_len, 0.0],
                          [0.0, rtrans_len]])
    
    rep["b"] = 0.5 * np.matrix(
        [[trans_len + rtrans_len, rtrans_len - trans_len],
         [rtrans_len - trans_len, trans_len + rtrans_len]]
    )

    return (rep, free2)

def psl2z_rep():
    # return the FSA for the presentation:
    # PSL(2,Z) = <a,b: a^2 = b^3 = 1>
    # and a representation into SL(2,R) known to be Anosov
    rep = Representation(["a", "b"])

    mu = 3.0
    D = np.matrix([[mu, 0], [0, 1/mu]])
    R2 = np.matrix([[0, -1], [1, 0]])
    R3 = np.matrix([[np.cos(np.pi / 3), -1 * np.sin(np.pi/3)],
                    [np.sin(np.pi / 3), np.cos(np.pi/3)]])
    
    rep["A"] = np.linalg.inv(D) * R3 * D
    rep["B"] = R2

    gdict = {1:{2:["a", "A"]}, 2:{1:["b"]}}
    return (rep, FSA(gdict))

def psl2z_multicones():
    #some precomputed multicones for the PSL(2,Z) representation given
    #above
    multicones = {}
    theta1 = np.pi/6
    theta2 = np.pi/8
    clusters = {
        1: np.array([[np.cos(theta1), np.sin(theta1)],
                     [np.cos(-1 * theta1), np.sin(-1 * theta1)]]),
        2: np.array([[np.cos(np.pi/2 + theta2), np.sin(np.pi/2 + theta2)],
                     [np.cos(np.pi/2 - theta2), np.sin(np.pi/2 - theta2)]])
    }

    return {v:Multicone(cluster)
            for v, cluster in clusters.iteritems()}

def genus2_rep():
    #return FSA for a genus 2 surface with presentation:
    #
    # <a,b,c,d | abABcdCD = 1>
    #
    # and the rep into SL(2,R) given by the regular octagon tiling of
    # H^2
    rep = Representation(["a", "b", "c", "d"])
    mats = np.load("g2s.npy")
    rep["a"] = so_to_psl(mats[0])
    rep["b"] = so_to_psl(mats[1])
    rep["c"] = so_to_psl(mats[2])
    rep["d"] = so_to_psl(mats[3])


    graph = FSA(load_fsa_dict("kbmag_files/genus2_surface.geowa"))
    graph.make_recurrent()

    return (rep, graph)

def hitchin_rep(dim):
    #same as above, but compose the rep into SL(2,R) with an irrep
    #into SL(n,R)
    rep, graph = genus2_rep()
    h_rep = Representation(["a", "b", "c", "d"])
    for g in "abcd":
        h_rep[g] = psl_irrep(rep[g], dim)
    return (h_rep, graph)
    
def transformed_cones(rep, graph, multicones):
    #for each vertex, return the images of the multicone associated to
    #that vertex under the transformations given by the outgoing edges
    cones = {}
    for v, cone in multicones.iteritems():
        for vt, w, label in graph.edges_out(v):
            mat = rep[rep.invert_gen(label)]
            cones[(vt,w,label)] = multicones[v].transform(mat)
    return cones

def get_angles(multicone):
    #given a multicone in R^2 (i.e. lifts of points in RP^1), get a
    #pair of bounding angles in S^1 for drawing purposes
    
    #assume that multicones take up <50% of RP^1
    points = sphere_lift(multicone.points)
    angles = np.arctan2(points[:,0], points[:,1]) * 180. / np.pi
    amin, amax = min(angles), max(angles)
    if amax - amin <= amin + 360. - amax:
        return (amin, amax)
    else:
        return (amax, amin + 360.)

def plot_multicones(multicones, radius = 1.0, ax = None, color='b'):
    #plot some multicones in RP^1
    unit_circle = plt.Circle((0,0), 1, fill=False)
    if ax is None:
        ax = plt.gca()
        ax.add_artist(unit_circle)
        ax.set_xlim(-1.2, 1.2)
        ax.set_ylim(-1.1, 1.2)
    
    for v, cone in multicones.items():
        theta1, theta2 = get_angles(cone)
        conewedge = Wedge((0,0), radius, theta1, theta2, color=color)
        ax.add_patch(conewedge)

    return ax
#        points = sphere_lift(cone.points).T
#        points = points / np.linalg.norm(points, axis = 0)
#        print points
#        plt.scatter(points[0], points[1], c=next(color))
    #plt.show()

def plot_clusters(clusters):
    color = iter(plt.cm.rainbow(np.linspace(0, 1, len(clusters))))
    for v, cluster in clusters.items():
        plt.scatter((cluster.T)[0], (cluster.T)[1], c=next(color))
    plt.show()
