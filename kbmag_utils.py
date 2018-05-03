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
                if vertex not in n_dict:
                    n_dict[vertex] = []
                n_dict[vertex].append(label)
                
        v_dict[i+1] = n_dict
    return v_dict

def run_kbmag(filename):
    autgroup = path.join(KBMAG_PATH, "autgroup")
    gpgeowa = path.join(KBMAG_PATH, "gpgeowa")


    subprocess.call([autgroup, filename])
    subprocess.call([gpgeowa, filename])

def load_fsa_dict(fsa_file):
    gap_dict = gap_record_parse.load_record_file(fsa_file)
    for fsa_dict in gap_dict.values():
        if fsa_dict["isFSA"] == "true":
            labels = fsa_dict["alphabet"]["names"]
            transitions = fsa_dict["table"]["transitions"]

            return build_dict(transitions, labels, to_filter = [0])
