import numpy as np
from src.pos_array import pos_array

def define_essential(mesh, problem, geometry, numbers, **kwargs):
    """
    Define essential degrees of freedom.

    Parameters:
        mesh: mesh structure
        problem: problem structure
        geometry: the 'geometry' to put essential degrees on
                  'nodes': in the nodes given by numbers
                  'points': in the points given by numbers
                  'curves': in the curves given by numbers
        numbers: an array of 'geometry' numbers
        **kwargs: Optional arguments
                  'physq': physical quantity number. default=1
                  'degfd': degree of freedom within the physical quantity. default=1
                  'iessp': index of previously defined essential degrees. Newly defined bc
                           will be added.

    Returns:
        iess: index of defined essential degrees
    """
    
    # Optional arguments
    physq = kwargs.get('physq', 1)
    degfd = kwargs.get('degfd', 1)
    iessp = kwargs.get('iessp', 0)
    add = iessp != 0

    # Define geometry
    if geometry == 'nodes':
        nodes = numbers
        nnodes = len(numbers)
    elif geometry == 'points':
        nodes = mesh.points[numbers]
        nnodes = len(numbers)
    elif geometry == 'curves':
        nnodes = sum([mesh.curves[curve].nnodes for curve in numbers])
        nodes = np.array([])
        for curve in numbers:
            np.append(nodes,mesh.curves[curve].nodes)
    else:
        raise ValueError(f"Invalid geometry: {geometry}")

    # Determine pos array
    pos = np.zeros(nnodes, dtype=int)
    ipos = 0

    for node in nodes:
        posn, ndof = pos_array(problem, node, physq=physq, order='ND')

        if degfd > ndof[0]:
            continue

        ipos += 1
        pos[ipos - 1] = posn[0][degfd - 1]  # Python 0-based indexing

    # Output iess
    if add:
        iess = np.unique(np.concatenate((pos[:ipos], [iessp])))
    else:
        iess = np.unique(pos[:ipos])

    return iess