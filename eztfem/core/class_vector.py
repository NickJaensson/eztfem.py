import numpy as np


# class definition for Vector objects
class Vector:

    # initialize attributes (will be filled elsewhere)
    def __init__(self, vec=0):
        self.vec = vec
        self.u = np.array([])

    # equivalence check for testing against Matlab code
    # NOTE: see the file NOTE_ON_COMPARING_ARRAYS.md for the use of np.squeeze
    def __eq__(self, other):
        check = [self.vec == other.vec,
                 np.allclose(np.squeeze(self.u), other.u, atol=1e-12, rtol=0)]
        return all(check)
