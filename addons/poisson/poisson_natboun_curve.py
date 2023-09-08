import numpy as np
from src.isoparametric_deformation_curve import isoparametric_deformation_curve

def poisson_natboun_curve(elem, coor, user, pos):
    """
    Boundary element for a natural boundary on a curve for the Poisson/diffusion equation:
    - alpha * dudn = h

    Parameters:
        elem (int): Element number.
        coor (np.ndarray): Coordinates of the nodes of the element.
                           Shape: (number of points in space, direction in space).
        user (dict): Dictionary containing user data.
            user.phi: Basis functions.
            user.dphi: Derivatives of basis functions.
            user.wg: Weights of the integration points.
            user.coorsys: Coordinate system (0 for Cartesian, 1 for axisymmetric).
            user.funcnr: Function number for the right-hand side h (flux).
            user.func: Function for the right-hand side.
                          Callable with signature: func(funcnr, x)
        pos (list): List of positions of degrees of freedom of each physq.

    Returns:
        elemvec (np.ndarray): The element vector.
    """

    # Set some values
    ninti = user.phi.shape[0]  # Number of integration points
    ndf = user.phi.shape[1]  # Number of degrees of freedom
    ndim = coor.shape[1]  # Dimension of space

    # Compute mapping of reference to real element
    dxdxi, curvel, normal = isoparametric_deformation_curve(coor, user.dphi)

    # Position of the integration points
    xg = user.phi.dot(coor)

    # Axisymmetric condition
    if user.coorsys == 1:
        curvel *= 2 * np.pi * xg[:, 1]

    # Compute element vector
    elemvec = np.zeros(ndf)

    if user.funcnr > 0:
        fg = np.zeros(ninti)
        for ip in range(ninti):
            fg[ip] = user.func(user.funcnr, xg[ip])
        elemvec = -user.phi.T.dot(fg * curvel * user.wg)

    return elemvec
