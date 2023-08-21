import numpy as np
from src.mesh_class import Mesh

def refcoor_nodal_points(mesh):
    """
    Return the reference coordinates of the nodal points.

    Parameters:
    - mesh (dict): Dictionary containing the mesh information. 
                   'elshape' key is required in the dictionary to determine the element shape.

    Returns:
    - x (np.ndarray): Reference coordinates of the nodal points.
    """

    elshape = mesh.elshape

    if elshape == 1:  # two-node line element
        x = np.array([-1, 1])
    elif elshape == 2:  # three-node line element
        x = np.array([-1, 0, 1])
    elif elshape == 3:  # three-node triangular element
        x = np.array([[0, 1, 0], [0, 0, 1]]).T
    elif elshape == 4:  # six-node triangular element
        x = np.array([[0, 0.5, 1, 0.5, 0, 0], [0, 0, 0, 0.5, 1, 0.5]]).T
    elif elshape == 5:  # four-node quadrilateral
        x = np.array([[-1, 1, 1, -1], [-1, -1, 1, 1]]).T
    elif elshape == 6:  # nine-node quadrilateral
        x = np.array([[-1, 0, 1, 1, 1, 0, -1, -1, 0], [-1, -1, -1, 0, 1, 1, 1, 0, 0]]).T
    elif elshape == 7:  # seven-node triangular element
        x = np.array([[0, 0.5, 1, 0.5, 0, 0, 1/3], [0, 0, 0, 0.5, 1, 0.5, 1/3]]).T
    elif elshape == 9:  # five-node quadrilateral
        x = np.array([[-1, 1, 1, -1, 0], [-1, -1, 1, 1, 0]]).T
    elif elshape == 10:  # four-node triangular element
        x = np.array([[0, 1, 0, 1/3], [0, 0, 1, 1/3]]).T
    else:
        raise ValueError(f"Invalid elshape = {elshape}")

    return x