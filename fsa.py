import copy

class FSA:
    def __init__(self, vert_dict):
        self._out_dict = copy.deepcopy(vert_dict)
        self._in_dict = self.in_dict()

    def in_dict(self):
        in_dict = {v:{} for v in self._out_dict}
        for v, out_neighbors in self._out_dict.iteritems():
            for w, label in out_neighbors.iteritems():
                in_dict[w][v] = label
        return in_dict

    def out_neighbors(self, vertex):
        return self._out_dict[vertex]

    def in_neighbors(self, vertex):
        return self._in_dict[vertex]

    def add_vertices(self, vertices):
        for v in vertices:
            if v not in self._out_dict:
                self._out_dict[v] = {}
                self._in_dict[v] = {}
                
    def add_edges(self, edges):
        for e in edges:
            tail, head, label = e
            self._out_dict[tail][head] = label
            self._in_dict[head][tail] = label

    def delete_vertices(self, vertices):
        for v in vertices:
            self.delete_vertex(v)

    def delete_vertex(self, vertex):
        for w in self.out_neighbors(vertex):
            self._in_dict[w].pop(vertex)
        for w in self.in_neighbors(vertex):
            self._out_dict[w].pop(vertex)
                
        self._out_dict.pop(vertex)
        self._in_dict.pop(vertex)
            
    def vertices(self):
        return self._out_dict.keys()

    def make_recurrent(self):
        still_pruning = True
        while still_pruning:
            still_pruning = False
            for v in self._out_dict.keys():
                if (len(self._out_dict[v]) == 0 or
                    len(self._in_dict[v]) == 0):
                    self.delete_vertex(v)
                    still_pruning = True
