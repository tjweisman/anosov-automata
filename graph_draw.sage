from collections import deque
from os import path
import subprocess


import gap_record_parse


KBMAG_PATH = "/home/teddy/math/tools/kbmag/bin"

def build_dict(v_list, labels, to_filter = []):
    v_dict = {}
    for i,neighbors in enumerate(v_list):
        n_dict = {}
        for label, vertex in zip(labels, neighbors):
            if vertex not in to_filter:
                n_dict[vertex] = label
                
        v_dict[i+1] = n_dict
    return v_dict

def run_kbmag(filename):
    autgroup = path.join(KBMAG_PATH, "autgroup")
    gpgeowa = path.join(KBMAG_PATH, "gpgeowa")

    print autgroup

    subprocess.call([autgroup, filename])
    subprocess.call([gpgeowa, filename])

def get_graph_from_fsa(fsa_file):
    gap_dict = gap_record_parse.load_record_file(fsa_file)
    for fsa_dict in gap_dict.values():
        if fsa_dict["isFSA"] == "true":
            labels = fsa_dict["alphabet"]["names"]
            transitions = fsa_dict["table"]["transitions"]

            vert_dict = build_dict(transitions, labels, to_filter = [0])
            return DiGraph(vert_dict, format = "dict_of_dicts")

def remove_long_paths(G, root, edge_ties=True):
    H = DiGraph([])
    H.add_vertices(G.vertices())

    #Dijkstra's algorithm!
    distance = {-1:False for v in G.vertices()}
    marked = {v:False for v in G.vertices()}

    marked[root] = True
    distance[root] = -1
    vertex_queue = deque([root])
    while len(vertex_queue) > 0:
        v = vertex_queue.popleft()
        
        to_visit = [w for w in G.neighbors_out(v) if not marked[w]]
        for w in to_visit:
            marked[w] = True
            distance[w] = distance[v] + 1

        short_nbrs = to_visit
        if edge_ties:
            short_nbrs = [w for w in G.neighbors_out(v)
                           if distance[w] == distance[v] + 1]
        
        H.add_edges([(v,w,G.edge_label(v,w)) for w in short_nbrs])
        vertex_queue.extend(to_visit)

    return H


