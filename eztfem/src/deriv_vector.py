import numpy as np
from .pos_array import pos_array
from .pos_array_vec import pos_array_vec
from .class_vector import Vector

def deriv_vector(mesh, problem, element, user, **kwargs):
    """
    Derive a vector of special structure.

    Parameters
    ----------
    mesh : Mesh
        Mesh structure.
    problem : Problem
        Problem structure.
    element : callable
        Function handle to the element function routine.
    user : Any
        User object to pass parameters and data to the element routine.
    **kwargs : dict, optional
        Optional arguments:
        vec : int, optional
            The vector number. Default is problem['nphysq'].
        order : str, optional
            The sequence order of the degrees of freedom on element level.
            'ND' : the most inner loop is over the degrees of freedom.
            'DN' : the most inner loop is over the nodal points.
            Default is 'DN'.
        posvectors : bool, optional
            Supply the position of vectors to the element routine. Default: False.

    Returns
    -------
    v : Vector
        Vector object containing the derived data.
    """

    # Set default optional arguments
    vec = kwargs.get('vec', problem.nphysq)  # first vec after nphysq (starting at zero)
    order = kwargs.get('order', 'DN')
    posvectors = kwargs.get('posvectors', 0)

    for kwarg in kwargs:
        if kwarg not in ['vec', 'order', 'posvectors']:
            raise ValueError(f'Invalid argument: {kwarg}')

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