import numpy as np
from ...core.isoparametric_deformation import isoparametric_deformation


def poisson_deriv(elem, coor, user, pos):
    """
    Compute the gradient of the velocity field for a given element.

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
        Gradient of the velocity field at all nodes, shape (n_dim * n_nodes,).

    Notes
    -----
    This function is intended for use with `deriv_vector` or `plot_sol`.

    """

    # Set some values
    ndim = coor.shape[1]
    nodalp = user.phi.shape[0]
    ndf = user.phi.shape[1]

    # Compute mapping of reference to real element
    F, Finv, detF = isoparametric_deformation(coor, user.dphi)

    # Compute derivative of the basis functions with respect to the real
    # coordinates
    dphidx = np.zeros((nodalp, ndf, ndim))
    for ip in range(nodalp):
        dphidx[ip, :, :] = user.dphi[ip, :, :].dot(Finv[ip, :, :])

    # Compute gradient vector
    gradu = np.zeros((nodalp, ndim))
    for j in range(ndim):
        gradu[:, j] = dphidx[:, :, j] @ user.u[pos[0]]

    elemvec = np.reshape(gradu, [ndim * nodalp], order='F')

    return elemvec
