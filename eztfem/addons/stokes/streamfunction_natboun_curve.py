import numpy as np
from ...core.isoparametric_deformation_curve import \
    isoparametric_deformation_curve


def streamfunction_natboun_curve(elem, coor, user, pos, posvec):
    """
    Compute the boundary element for a natural boundary on a curve for the
    streamfunction equation (Poisson equation): dpsidn = -v * nx + u * ny
    since dpsidx = -v, dpsidy = u.

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
    posvec : list of ndarray
        Positions of the degrees of freedom of each vector.

    Returns
    -------
    elemmat : ndarray
        The element matrix.
    elemvec : ndarray
        The element vector.

    Notes
    -----
    This function must be called in `build_system` using `posvectors=1`.

    """

    # Set some values
    ninti = user.phi.shape[0]  # Number of integration points
    ndf = user.phi.shape[1]  # Number of degrees of freedom
    ndim = coor.shape[1]  # Dimension of space

    # Compute mapping of reference to real element
    dxdxi, curvel, normal = isoparametric_deformation_curve(coor, user.dphi)

    # Position of the integration points
    xg = user.phi @ coor

    if user.coorsys == 1:
        # Axisymmetric
        curvel = 2 * np.pi * xg[:, 1] * curvel

    # Compute element vector
    u = user.v[posvec[1]]  # Velocity
    tmp = u.reshape((ndf, ndim), order='F')

    gradpsi = np.zeros((ninti, 2))

    gradpsi[:, 0] = -np.dot(user.phi, tmp[:, 1])
    gradpsi[:, 1] = np.dot(user.phi, tmp[:, 0])

    work = np.zeros(ninti)

    for ip in range(ninti):
        work[ip] = 0
        for i in range(2):
            work[ip] += normal[ip, i] * gradpsi[ip, i]

    elemvec = np.dot(user.phi.T, work * curvel * user.wg)

    return elemvec
