import numpy as np
from .pos_array import pos_array


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
    **kwargs : optional
        Optional arguments:
        - funcnr : int, optional
            Function number for func, i.e. func(funcnr, x). Default is 1.
        - physq : int, optional
            Physical quantity number. Default is 1.
        - degfd : int, optional
            Degree of freedom within the physical quantity. Default is 1.
        - f : numpy.ndarray, optional
            An existing system vector.

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
