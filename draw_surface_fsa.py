import sys
import matplotlib.pyplot as plt
from math import sin, cos, pi

from kbmag_utils import load_fsa_dict
import fsa

def arrow(tail, head, color="k"):
    ax = plt.axes()
    ax.arrow(tail[0], tail[1],
             head[0] - tail[0],
             head[1] - tail[1],
             head_width = 0.1,
             head_length = 0.1,
             head_starts_at_zero = True,
             fc = color,
             ec = color)


def main(args):
    G = fsa.FSA(load_fsa_dict("kbmag_files/genus2_surface.geowa"))
    H = fsa.remove_long_paths(G, 1)
    levels = {0:[1],
              1:[4,6,9,7,8,2,5,3],
              2:[15,14,18,19,24,25,21,20,23,22,10,11,17,16,12,13],
              3:[30,34,35,40,41,37,36,39,38,26,27,33,32,28,29,31],
              4:[46,47,49,48,42,43,44,45]}

    positions = {}
    delta = {0:0,
             1:0,
             2:-1 * pi / 16,
             3: pi / 16,
             4: pi / 8}

    plt.axis("off")
    for level, vertices in levels.iteritems():
        n = len(vertices)
        for i, v in enumerate(vertices):
            positions[v] = (float(level) * cos(i * 2 * pi / n + delta[level]),
                            float(level) * sin(i * 2 * pi / n + delta[level]))
            
        x = [float(level) * cos(i * 2 * pi / n + delta[level]) for i in range(n)]
        y = [float(level) * sin(i * 2 * pi / n + delta[level]) for i in range(n)]

        plt.plot(x,y, "o")

    for v,w in G.edges():
        if w not in H.neighbors_out(v):
            arrow(positions[v], positions[w], "gray")
        
    for v,w in H.edges():
        arrow(positions[v], positions[w], "black")

    

    plt.show()
    
    

if __name__ == "__main__":
    main(sys.argv[1:])
