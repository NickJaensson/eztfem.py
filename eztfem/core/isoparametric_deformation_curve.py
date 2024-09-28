import numpy as np


def isoparametric_deformation_curve(x, dphi):
    """
    Isoparametric deformation of curved line elements.

    Parameters
    ----------
    x : numpy.ndarray
        Coordinates of the nodes, representing the position of the unknowns for
        phi and dphi.
        Shape: (number of points in space, direction in space).
    dphi : numpy.ndarray
        Derivative of the shape function with respect to the reference
        coordinate.
        Shape: (number of points in space, number of unknowns, 1).

    Returns
    -------
    dxdxi : numpy.ndarray
        Derivative of the coordinates x with respect to the reference
        coordinate (tangential vector).
        Shape: (number of points in space, direction in space).
    curvel : numpy.ndarray
        The length of the vector dxdxi.
        Shape: (number of points in space).
    normal : numpy.ndarray
        Unit normal vector to the curve.
        Shape: (number of points in space, direction in space).

    Examples
    --------
    >>> dxdxi, curvel, normal = isoparametric_deformation_curved_line(x, dphi)
    """

    npts, ndim = dphi.shape[0], x.shape[1]

    dxdxi = np.zeros((npts, ndim))
    curvel2 = np.zeros(npts)

    # Compute dxdxi and curvel2
    for j in range(ndim):
        dxdxi[:, j] = dphi[:, :, 0] @ x[:, j]
        curvel2 += dxdxi[:, j]**2

    curvel = np.sqrt(curvel2)

    normal = np.zeros((npts, 2))

    if ndim == 2:
        normal[:, 0] = dxdxi[:, 1] / curvel
        normal[:, 1] = -dxdxi[:, 0] / curvel

    return dxdxi, curvel, normal
