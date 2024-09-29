import numpy as np
from ...core.isoparametric_deformation import isoparametric_deformation


def stokes_deriv(elem, coor, user, pos):
    """
    Compute the derivative of the velocity field for post-processing.

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
        The computed derivative at all nodes.

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

    # Get velocity vector
    u = user.u[pos[0]]
    uvector = u.reshape(ndf, ndim, order='F')  # use Fortran numbering

    # Compute velocity gradient tensor
    gradu = np.zeros((nodalp, ndim, ndim))

    for j in range(ndim):
        gradu[:, :, j] = dphidx[:, :, j].dot(uvector)

    # Element vector
    if user.comp == 0:
        elemvec = gradu[:, 0, 0]
    elif user.comp == 1:
        elemvec = gradu[:, 0, 1]
    elif user.comp == 2:
        elemvec = gradu[:, 1, 0]
    elif user.comp == 3:
        elemvec = gradu[:, 1, 1]
    elif user.comp == 4:
        elemvec = gradu[:, 1, 0] - gradu[:, 0, 1]
    elif user.comp == 5:
        if user.coorsys == 0:
            raise ValueError('gradu in theta direction not applicable for \
                             coorsys=0')
        smallr = coor[:, 1] < 1e-10
        elemvec = np.zeros(nodalp)
        elemvec[smallr] = gradu[smallr, 1, 1]
        elemvec[~smallr] = uvector[~smallr, 1] / coor[~smallr, 1]
    elif user.comp == 6:
        elemvec = gradu[:, 0, 0] + gradu[:, 1, 1]
        if user.coorsys == 1:
            smallr = coor[:, 1] < 1e-10
            elemvec[smallr] += gradu[smallr, 1, 1]
            elemvec[~smallr] += uvector[~smallr, 1] / coor[~smallr, 1]
    elif user.comp == 7:
        elemvec = np.zeros(nodalp)
        for ip in range(nodalp):
            gu = gradu[ip, :, :]
            elemvec[ip] = np.sum((gu + gu.T) ** 2) / 2
            if user.coorsys == 1:
                if coor[ip, 1] < 1e-10:
                    elemvec[ip] += 2 * (gradu[ip, 1, 1]) ** 2
                else:
                    elemvec[ip] += 2 * (uvector[ip, 1] / coor[ip, 1]) ** 2
        elemvec = np.sqrt(elemvec)
    else:
        raise ValueError(f"Invalid user.comp: {user.comp}")

    return elemvec
