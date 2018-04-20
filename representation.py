import numpy as np
import re

class RepresentationException:
    pass

class Representation:
    #technically, this is just a free group representation

    @staticmethod
    def invert_gen(generator):
        if re.match("[a-z]", generator):
            return generator.upper()
        else:
            return generator.lower()

    @property
    def dim(self):
        return self._dim
    
    def __init__(self, generator_names=[]):
        self.generators = {name[0].lower():None
                           for name in generator_names}
        
        for gen, image in self.generators.items():
            self.generators[gen.upper()] = None
            
        self._dim = None

    def __getitem__(self, generator):
        if len(generator) == 1:
            return self.generators[generator]
        elif len(generator) > 1:
            return self[generator[0]] * self[generator[1:]]
        else:
            return np.identity(self._dim)

    def __setitem__(self, generator, matrix):
        shape = matrix.shape

        if self._dim is None:
            self._dim = shape[0]
        if shape[0] != shape[1]:
            raise RepresentationException("use square matrices")
        if shape[0] != self._dim:
            raise RepresentationException("use matrices of matching dimensions")
        
        self.generators[generator] = matrix
        self.generators[self.invert_gen(generator)] = np.linalg.inv(matrix)
