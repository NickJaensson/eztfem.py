import numpy as np
from .pos_array import pos_array, pos_array_vec


class Vector:
    """
    A class used to represent a Vector.

    Attributes
    ----------
    vec : int
        The vector value.
    u : np.ndarray
        An array to store vector data.

    """

    def __init__(self, vec=0):
        """
        Initializes the Vector object with the given attributes. The stored
        vecor data (u) is initialized as an empty numpy array (np.array([])).

        Parameters for initialization
        -----------------------------
        vec : int or float, optional
            The vector value (default is 0).

        """
        self.vec = vec
        self.u = np.array([])

    def __eq__(self, other):
        """
        Checks equivalence of two Vector objects (overloads == sign).

        Parameters
        ----------
        other : Vector
            The other Vector object to compare with.

        Returns
        -------
        bool
            True if the Vector objects are equivalent, False otherwise.

        Notes
        -----
        See NOTE_ON_COMPARING_ARRAYS.md for the use of np.squeeze.

        """
        check = [self.vec == other.vec,
                 np.allclose(np.squeeze(self.u), other.u, atol=1e-12, rtol=0)]
        return all(check)


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

    Keyword arguments
    -----------------
    vec : int, optional
        The vector number. Default is problem['nphysq'].
    order : str, optional
        The sequence order of the degrees of freedom on element level.
        'ND' : the most inner loop is over the degrees of freedom.
        'DN' : the most inner loop is over the nodal points.
        Default is 'DN'.
    posvectors : bool, optional
        Supply the position of vectors to the element routine.
        Default: False.

    Returns
    -------
    v : Vector
        Vector object containing the derived data.

    """

    # Set default optional arguments
    # NOTE: default vec is first vec after nphysq (starting at zero)
    vec = kwargs.get('vec', problem.nphysq)
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
        posvec, _ = pos_array_vec(problem, mesh.topology[:, elem].T, vec=vec,
                                  order=order)
        posv = np.array(posvec).flatten()
        coor = mesh.coor[mesh.topology[:, elem], :]

        if posvectors:
            posvec, _ = pos_array_vec(problem, mesh.topology[:, elem].T,
                                      order=order)
            elemvec = element(elem, coor, user, pos, posvec)
        else:
            elemvec = element(elem, coor, user, pos)

        v.u[posv] += elemvec
        w[posv] += 1.0

    # Apply averaging by the total weight
    in_ = w > 0
    v.u[in_] = v.u[in_] / w[in_]

    return v


def fill_system_vector(mesh, problem, geometry, numbers, func, **kwargs):
    """
    Fill system vector.

    NOTE: If f is present, this function modifies the input arguments `f`
          in place. If f is not present, f is returned by the function.

    Parameters
    ----------
    mesh : object
        Mesh structure.
    problem : object
        Problem structure.
    geometry : str
        The 'geometry' to fill degrees on:
        - 'nodes': in the nodes given by numbers
        - 'points': in the points given by numbers
        - 'curves': in the curves given by numbers
    numbers : array_like
        An array of 'geometry' numbers.
    func : callable
        Scalar function for filling. Only argument is x, a coordinate vector.

    Keyword arguments
    -----------------
    funcnr : int, optional
        Function number for func, i.e. func(funcnr, x). Default is 1.
    physq : int, optional
        Physical quantity number. Default is 1.
    degfd : int, optional
        Degree of freedom within the physical quantity. Default is 1.
    f : numpy.ndarray, optional
        An existing system vector which will be modified in-place.

    Returns
    -------
    f : numpy.ndarray (optional)
        Filled system vector. Will only be returned if f is not present in the
        function call (if it is present, f is modified in place).

    """

    funcnr = kwargs.get('funcnr', 0)
    physq = kwargs.get('physq', 0)
    degfd = kwargs.get('degfd', 0)
    f = kwargs.get('f', None)

    f_present = f is not None

    if not f_present:
        f = np.zeros(problem.numdegfd)

    # Convert numbers to a list if an int is supplied
    if isinstance(numbers, (int, np.integer)):
        numbers = [numbers]

    # Define geometry
    if geometry == 'nodes':
        nodes = numbers
    elif geometry == 'points':
        nodes = mesh.points[numbers]
    elif geometry == 'curves':
        # nnodes = sum([mesh.curves[curve].nnodes for curve in numbers])
        nodes = np.array([], dtype=int)
        for curve in numbers:
            nodes = np.append(nodes, mesh.curves[curve].nodes)
    else:
        raise ValueError(f"Invalid geometry: {geometry}")

    # Fill vector
    for node in nodes:
        posn, ndof = pos_array(problem, node, physq=physq, order='ND')

        if degfd > ndof[0]:
            continue

        f[posn[0][degfd]] = func(funcnr, mesh.coor[node])

    if not f_present:
        return f
