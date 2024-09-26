import numpy as np
from .pos_array import pos_array

def define_essential(mesh, problem, geometry, numbers, **kwargs):
    """
    Get the indices of the essential degrees of freedom.

    Parameters
    ----------
    mesh : Mesh
        Mesh object.
    problem : Problem
        Problem object.
    geometry : str
        The 'geometry' to put essential degrees on:
        'nodes' : in the nodes given by numbers.
        'points' : in the points given by numbers.
        'curves' : in the curves given by numbers.
    numbers : array_like
        An array of 'geometry' numbers.
    **kwargs : dict, optional
        Optional arguments:
        physq : int, optional
            Physical quantity number. Default is 0.
        degfd : int, optional
            Degree of freedom within the physical quantity. Default is 0.
        iessp : int, optional
            Index of previously defined essential degrees. Newly defined boundary 
            conditions will be added. Default is 0.

    Returns
    -------
    iess : numpy.ndarray
        Index of defined essential degrees.
    """
    
    # Optional arguments
    physq = kwargs.get('physq', 0)
    degfd = kwargs.get('degfd', 0)
    iessp = kwargs.get('iessp', 0)
    add = 'iessp' in kwargs

    # Convert numbers to a list if an int is supplied
    if isinstance(numbers, (int, np.integer)):
        numbers = [numbers]

    # Define geometry
    if geometry == 'nodes':
        nodes = numbers
        nnodes = len(numbers)
    elif geometry == 'points':
        nodes = mesh.points[numbers]
        nnodes = len(numbers)
    elif geometry == 'curves':
        nnodes = sum([mesh.curves[curve].nnodes for curve in numbers])
        nodes = np.array([],dtype=int)
        for curve in numbers:
            nodes = np.append(nodes,mesh.curves[curve].nodes)
    else:
        raise ValueError(f"Invalid geometry: {geometry}")

    # Determine pos array
    pos = np.zeros(nnodes, dtype=int)
    ipos = 0

    for node in nodes:
        posn, ndof = pos_array(problem, node, physq=physq, order='ND')

        if degfd > ndof[0]:
            continue

        pos[ipos] = posn[0][degfd]  # Python 0-based indexing
        ipos += 1

    # Output iess
    if add:
        # Convert iessp to a numpy array if an int is supplied
        if isinstance(iessp, (int, np.integer)):
            iessp = np.array([iessp])
        iess = np.unique(np.concatenate((pos[:ipos], iessp)))
    else:
        iess = np.unique(pos[:ipos])

    return iess