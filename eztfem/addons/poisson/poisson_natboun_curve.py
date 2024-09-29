import numpy as np
from ...core.isoparametric_deformation_curve import \
    isoparametric_deformation_curve


def poisson_natboun_curve(elem, coor, user, pos):
    """
    Boundary element for a natural boundary on a curve for the
    Poisson/diffusion equation: alpha * dudn = h

    Parameters
    ----------
    elem : int
        Element number.
    coor : ndarray
        Coordinates of the nodes of the element, shape (n_points, n_dim).
    user : User
        User object containing shape function, Gauss points, parameters
    pos : list of ndarray
        Positions of the degrees of freedom of each physical quantity.

    Returns
    -------
    elemvec : ndarray
        The element vector.

    Notes
    -----
    This function computes the element vector for a natural boundary on a curve
    for the Poisson/diffusion equation using isoparametric finite element
    methods.

    """

    # Set some values
    ninti = user.phi.shape[0]  # Number of integration points
    ndf = user.phi.shape[1]    # Number of degrees of freedom

    # Compute mapping of reference to real element
    dxdxi, curvel, normal = isoparametric_deformation_curve(coor, user.dphi)

    # Position of the integration points
    xg = user.phi @ coor

    # Axisymmetric condition
    if user.coorsys == 1:
        curvel *= 2 * np.pi * xg[:, 1]

    # Compute element vector
    elemvec = np.zeros(ndf)

    if user.funcnr > 0:
        fg = np.zeros(ninti)
        for ip in range(ninti):
            fg[ip] = user.func(user.funcnr, xg[ip])
        elemvec = -user.phi.T @ (fg * curvel * user.wg)

    return elemvec
