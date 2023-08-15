import numpy as np
from src.pos_array import pos_array

def fill_system_vector(mesh, problem, geometry, numbers, func, **kwargs):
    """
    Fill system vector.

    Parameters:
        mesh: mesh structure
        problem: problem structure
        geometry: the 'geometry' to fill degrees on:
                  'nodes': in the nodes given by numbers
                  'points': in the points given by numbers
                  'curves': in the curves given by numbers
        numbers: an array of 'geometry' numbers
        func: scalar function for filling. Only argument is x, a coordinate vector.
        **kwargs: Optional arguments
                  'funcnr': function number for func, i.e. func(funcnr,x). default=1
                  'physq': physical quantity number. default=1
                  'degfd': degree of freedom within the physical quantity. default=1
                  'fin': an existing system vector

    Returns:
        f: filled system vector. If fin is present, f is a copy of fin with
           newly filled values where specified. Otherwise the non-filled part
           is initialized with zero.
    """
    
    funcnr = kwargs.get('funcnr', 0)
    physq = kwargs.get('physq', 0)
    degfd = kwargs.get('degfd', 0)
    fin = kwargs.get('fin', None)

    if fin is None:
        f = np.zeros(problem.numdegfd)
    else:
        f = np.array(fin) # copy the array

    # Convert numbers to a list if an int is supplied
    if isinstance(numbers, (int, np.integer)):
        numbers = [numbers]

    # Define geometry
    if geometry == 'nodes':
        nodes = numbers
    elif geometry == 'points':
        nodes = mesh.points[numbers]
    elif geometry == 'curves':
        nnodes = sum([mesh.curves[curve].nnodes for curve in numbers])
        nodes = np.array([],dtype=int)
        for curve in numbers:
            nodes = np.append(nodes,mesh.curves[curve].nodes)
    else:
        raise ValueError(f"Invalid geometry: {geometry}")

    # Fill vector
    for node in nodes:
        posn, ndof = pos_array(problem, node, physq=physq, order='ND')
        
        if degfd > ndof[0]:
            continue

        f[posn[0][degfd]] = func(funcnr, mesh.coor[node])

    return f
