'''
Element routines for the streamfunction equation.
'''
import numpy as np
from ...core.shapefunc import isoparametric_deformation, \
    isoparametric_deformation_curve


def streamfunction_elem(elem, coor, user, pos, posvec):
    """
    Compute the element matrix and vector for the streamfunction equation:
    - nabla^2 psi = omega, where omega is derived from the velocity vector.

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
    ninti = user.phi.shape[0]  # number of integration points
    ndf = user.phi.shape[1]  # number of degrees of freedom
    ndim = coor.shape[1]  # dimension of space

    # compute mapping of reference to real element
    _, fmat_inv, det_fmat = isoparametric_deformation(coor, user.dphi)

    # Position of the integration points
    xg = user.phi @ coor

    if user.coorsys == 1:
        # Axisymmetric
        det_fmat = 2 * np.pi * xg[:, 1] * det_fmat

    # compute derivative of the basis functions with respect to the real
    # coordinates
    dphidx = np.zeros((ninti, ndf, ndim))

    for ip in range(ninti):
        dphidx[ip, :, :] = user.dphi[ip, :, :] @ fmat_inv[ip, :, :]

    # Compute element matrix
    elemmat = np.zeros((ndf, ndf))
    work = np.zeros(ninti)

    for i in range(ndf):
        for j in range(i, ndf):
            work[:] = 0
            for ip in range(ninti):
                for k in range(ndim):
                    work[ip] += dphidx[ip, i, k] * dphidx[ip, j, k]
            elemmat[i, j] = np.sum(work * det_fmat * user.wg)
            elemmat[j, i] = elemmat[i, j]  # Symmetry

    # Compute element vector
    u = user.v[posvec[1]]  # Velocity
    tmp = u.reshape((ndf, ndim), order='F')

    dudy = np.dot(dphidx[:, :, 1], tmp[:, 0])
    dvdx = np.dot(dphidx[:, :, 0], tmp[:, 1])

    if user.coorsys == 1:
        work = np.dot(user.phi, u[:ndf])
        elemvec = 2 * np.pi * np.dot(user.phi.T, (xg[:, 1] * (dvdx - dudy)
                                                  - work) * det_fmat * user.wg)
    else:
        elemvec = np.dot(user.phi.T, (dvdx - dudy) * det_fmat * user.wg)

    return elemmat, elemvec


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
    _, curvel, normal = isoparametric_deformation_curve(coor, user.dphi)

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
