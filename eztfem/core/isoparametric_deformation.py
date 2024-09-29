import numpy as np


def isoparametric_deformation(x, dphi):
    """
    Isoparametric deformation of an element.

    Parameters
    ----------
    x : numpy.ndarray
        Coordinates of the nodes, representing the position of the unknowns
        for phi and dphi. x[i, j] with i being the point in space and j the
        direction in space.
    dphi : numpy.ndarray
        Derivatives of the shape function with respect to the reference
        coordinates.
        dphi[i, j, k], with i being the point in space, j the unknown, and k
        the direction in space.

    Returns
    -------
    F : numpy.ndarray
        Deformation gradient matrix between the reference element and the
        actual element.
        F[i, j, m] means that at point i, the transformation is:

        .. math::
            F[j, m] = \\frac{d x_j}{d \\xi_m}
    Finv : numpy.ndarray
        The inverse of F.
    detF : numpy.ndarray
        The determinant of F (Jacobian). detF[i] gives the determinant at
        point i.

    Examples
    --------
    >>> F, Finv, detF = isoparametric_deformation(x, dphi)

    """

    npts, ndim = dphi.shape[0], x.shape[1]

    F = np.zeros((npts, ndim, ndim))
    Finv = np.zeros((npts, ndim, ndim))
    detF = np.zeros(npts)

    # compute F
    for j in range(ndim):
        F[:, :, j] = dphi[:, :, j] @ x

    # compute detF and Finv
    if ndim == 1:
        # 1D
        for ip in range(npts):
            detF[ip] = F[ip, 0, 0]
            Finv[ip, 0, 0] = 1 / F[ip, 0, 0]
    elif ndim == 2:
        # 2D
        for ip in range(npts):
            detF[ip] = F[ip, 0, 0] * F[ip, 1, 1] - F[ip, 0, 1] * F[ip, 1, 0]
            Finv[ip, 0, 0] = F[ip, 1, 1] / detF[ip]
            Finv[ip, 0, 1] = -F[ip, 0, 1] / detF[ip]
            Finv[ip, 1, 0] = -F[ip, 1, 0] / detF[ip]
            Finv[ip, 1, 1] = F[ip, 0, 0] / detF[ip]

    return F, Finv, detF
