import numpy as np
from eztfem.src.pos_array import pos_array
from eztfem.src.pos_array_vec import pos_array_vec

def integrate_boundary_elements(mesh, problem, element, user, **kwargs):
    # Optional arguments
    curve = kwargs.get('curve', 0)
    order = kwargs.get('order', 'DN')
    posvectors = kwargs.get('posvectors', 0)

    if curve == 0:
        raise ValueError('Argument curve is missing.')

    # Start summing loop over elements
    resultsum = 0

    for elem in range(mesh.curves[curve].nelem ):
        # Positions in global system
        pos, _ = pos_array(problem, mesh.curves[curve].topology[:, elem, 1].T, order=order)
        coor = mesh.coor[mesh.curves[curve].topology[:, elem, 1], :]

        if posvectors:
            posvec, _ = pos_array_vec(problem, mesh.curves[curve].topology[:, elem, 1].T, order=order)
            resultsum += element(elem, coor, user, pos, posvec)
        else:
            resultsum += element(elem, coor, user, pos)

    return resultsum