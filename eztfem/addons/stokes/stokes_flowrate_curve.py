import numpy as np
from ...core.isoparametric_deformation_curve import \
      isoparametric_deformation_curve


def stokes_flowrate_curve(elem, coor, user, pos):
    """
    Compute the flowrate through a curve for boundary elements.

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
    flowrate : float
        The computed element flowrate.

    Notes
    -----
    This function must be called in `integrate_boundary_curve` using 
    `posvectors=0` (default).

    """

    # Function implementation goes here
    pass

    # Set some values
    ninti = user.phi.shape[0]  # Number of integration points
    ndf = user.phi.shape[1]    # Number of velocity dofs per spatial direction
    ndim = coor.shape[1]  # Dimension of space

    # Compute mapping of reference to real element
    dxdxi, curvel, normal = isoparametric_deformation_curve(coor, user.dphi)

    # Position of the integration points
    xg = np.dot(user.phi, coor)

    if user.coorsys == 1:
        # Axisymmetric
        curvel = 2 * np.pi * xg[:, 1] * curvel

    # Get velocity vector
    u = user.u[pos[0]]  # Velocity
    uvector_g = np.dot(user.phi, np.reshape(u, (ndf, ndim), order='F'))

    # Normal velocity
    un = np.zeros(ninti)

    for ip in range(ninti):
        un[ip] = np.dot(uvector_g[ip, :], normal[ip, :])

    # Compute flow rate for this element
    flowrate = np.sum(un * curvel * user.wg)

    return flowrate
