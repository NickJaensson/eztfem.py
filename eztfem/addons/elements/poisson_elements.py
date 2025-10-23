import numpy as np

from ...core.shapefunc import isoparametric_deformation, \
    isoparametric_deformation_curve
from .tp import NatbounCurveUser, PoissonElemUser, VelocityAwareUser



# FIXME: Docstring suggests pos is a list[ndarray], internal typing all the way
#        down to pos_array.py suggests these should take list[list[int]].
def poisson_elem(elem: int, coor: np.ndarray, user: PoissonElemUser, pos: list[list[int]]):
    """
    Element routine for the Poisson equation: - alpha nabla^2 u = f

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
    elemmat : ndarray
        Element matrix, shape (n_df, n_df).
    elemvec : ndarray
        Element vector, shape (n_df,).

    Notes
    -----
    This function computes the element matrix and vector for the Poisson
    equation using isoparametric finite element methods.

    """

    ninti = user.phi.shape[0]  # number of integration points
    ndf = user.phi.shape[1]    # number of degrees of freedom
    ndim = coor.shape[1]       # dimension of space

    # Compute mapping of reference to real element
    _F, Finv, detF = isoparametric_deformation(coor, user.dphi)

    # Position of the integration points
    xg = np.dot(user.phi, coor)

    if user.coorsys == 1:
        # axisymmetric
        detF = 2 * np.pi * xg[:, 1] * detF

    # Compute derivative of the basis functions wrt real coordinates
    dphidx = np.zeros((ninti, ndf, ndim))
    ndr = Finv.shape[1]

    for ip in range(ninti):
        dphidx[ip] = np.reshape(user.dphi[ip],
                                (ndf, ndr)) @ np.reshape(Finv[ip], (ndr, ndim))

    # Compute element matrix
    elemmat = np.zeros((ndf, ndf))
    work = np.zeros(ninti)

    for i in range(ndf):
        for j in range(i, ndf):
            work[:] = 0
            for ip in range(ninti):
                for k in range(ndim):
                    work[ip] += dphidx[ip, i, k] * dphidx[ip, j, k]
            elemmat[i, j] = user.alpha * np.sum(work * detF * user.wg)
            elemmat[j, i] = elemmat[i, j]  # symmetry

    # Compute element vector
    elemvec = np.zeros(ndf)

    if user.funcnr > 0:
        fg = np.zeros(ninti)

        for ip in range(ninti):
            fg[ip] = user.func(user.funcnr, xg[ip])

        for i in range(ndf):
            elemvec[i] = np.sum(fg * user.phi[:, i] * detF * user.wg)

    return elemmat, elemvec


def poisson_deriv(elem: int, coor: np.ndarray, user: VelocityAwareUser,
                  pos: list[list[int]]):
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
    print(pos)

    # Set some values
    ndim = coor.shape[1]
    nodalp = user.phi.shape[0]
    ndf = user.phi.shape[1]

    # Compute mapping of reference to real element
    _F, Finv, _detF = isoparametric_deformation(coor, user.dphi)

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


def poisson_natboun_curve(elem: int, coor: np.ndarray, user: NatbounCurveUser, pos: list[list[int]]):
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
    _dxdxi, curvel, _normal = isoparametric_deformation_curve(coor, user.dphi)

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
