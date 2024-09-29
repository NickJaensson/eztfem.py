import numpy as np


def refcoor_nodal_points(mesh):
    """
    Find the reference coordinates of the nodal points.

    Parameters
    ----------
    mesh : Mesh
        Mesh object.

    Returns
    -------
    x : numpy.ndarray
        Reference coordinates of the nodal points.

    """

    elshape = mesh.elshape

    match elshape:
        case 1:  # two-node line element
            x = np.array([-1, 1])
        case 2:  # three-node line element
            x = np.array([-1, 0, 1])
        case 3:  # three-node triangular element
            x = np.array([[0, 1, 0], [0, 0, 1]]).T
        case 4:  # six-node triangular element
            x = np.array([[0, 0.5, 1, 0.5, 0, 0], [0, 0, 0, 0.5, 1, 0.5]]).T
        case 5:  # four-node quadrilateral
            x = np.array([[-1, 1, 1, -1], [-1, -1, 1, 1]]).T
        case 6:  # nine-node quadrilateral
            x = np.array([[-1, 0, 1, 1, 1, 0, -1, -1, 0],
                          [-1, -1, -1, 0, 1, 1, 1, 0, 0]]).T
        case 7:  # seven-node triangular element
            x = np.array([[0, 0.5, 1, 0.5, 0, 0, 1/3],
                          [0, 0, 0, 0.5, 1, 0.5, 1/3]]).T
        case 9:  # five-node quadrilateral
            x = np.array([[-1, 1, 1, -1, 0], [-1, -1, 1, 1, 0]]).T
        case 10:  # four-node triangular element
            x = np.array([[0, 1, 0, 1/3], [0, 0, 1, 1/3]]).T
        case _:
            raise ValueError(f"Invalid elshape = {elshape}")

    return x
