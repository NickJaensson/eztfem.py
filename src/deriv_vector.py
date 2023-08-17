import numpy as np
from src.pos_array import pos_array
from src.pos_array_vec import pos_array_vec

# class definition for Mesh objects
class Vector:

    # initialize attributes (will be filled elsewhere)
    def __init__(self,vec=0):
        self.vec = vec
        self.u = np.array([])
        
    # equivalence check for testing against Matlab code
    # NOTE: see the matlab file create_pytfem_tests.m for the use of np.squeeze
    def __eq__(self, other):
       check = [self.vec == other.vec,
                np.allclose(np.squeeze(self.u),other.u,atol=1e-12,rtol=0)]
       return all(check)

def deriv_vector(mesh, problem, element, user, **kwargs):
    """
    Derive a vector of special structure.

    Parameters:
    - mesh (dict): Mesh structure.
    - problem (dict): Problem structure.
    - element (func): Function handle to the element function routine.
    - user: Can be used by the user for transferring data to the element routine.

    Keyword arguments:
    - vec (int): The vector number. Default is problem['nphysq'] + 1.
    - order (str): The sequence order of the degrees of freedom on element level.
                   'ND': the most inner loop is over the degrees of freedom.
                   'DN': the most inner loop is over the nodal points.
                   Default is 'DN'.
    - posvectors (int): Supply the position of vectors to the element routine. Default is 0.

    Returns:
    - v (dict): The vector. Contains 'vec' and 'u'.
    """

    # Default values for optional arguments
    vec = problem.nphysq # first vec after nphysq (starting at zero)
    order = 'DN'
    posvectors = 0

    # Process optional arguments
    for key, value in kwargs.items():
        if key == 'vec':
            vec = value
        elif key == 'order':
            order = value
        elif key == 'posvectors':
            posvectors = value
        else:
            raise ValueError(f"Invalid option: {key}")

    v = Vector(vec=vec)
    n = problem.vec_numdegfd[vec]
    v.u = np.zeros(n)
    w = np.zeros(n)  # weights

    # Start assemble loop over elements
    for elem in range(mesh.nelem):
        # Positions in the global system
        pos, _ = pos_array(problem, mesh.topology[:, elem].T, order=order)
        posvec, _ = pos_array_vec(problem, mesh.topology[:, elem].T, vec=vec, order=order)
        posv = np.array(posvec).flatten()
        coor = mesh.coor[mesh.topology[:,elem],:]

        if posvectors:
            posvec, _ = pos_array_vec(problem, mesh.topology[:, elem].T, order=order)
            elemvec = element(elem, coor, user, pos, posvec)
        else:
            elemvec = element(elem, coor, user, pos)

        v.u[posv] += elemvec
        w[posv] += 1.0

    # Apply averaging by the total weight
    in_ = w > 0
    v.u[in_] = v.u[in_] / w[in_]

    return v