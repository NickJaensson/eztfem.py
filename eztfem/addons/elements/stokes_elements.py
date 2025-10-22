import typing

import numpy as np
from ...core.shapefunc import isoparametric_deformation, \
    isoparametric_deformation_curve

if typing.TYPE_CHECKING:
    from ...core.user import User


# FIXME: Docstring suggests pos is a list[ndarray], internal typing all the way
#        down to pos_array.py suggests these should take list[list[int]].
def stokes_elem(elem: int, coor: np.ndarray, user: "User", pos: list[list[int]]):
    """
    Element routine for the Poisson equation:
    - nabla.( mu (nabla u+nabla u^T) ) + nabla p = f and nabla.u = 0

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
    This function computes the element matrix and vector for the Stokes
    equation using isoparametric finite element methods.

    """

    ndim = coor.shape[1]  # dimension of space
    ninti = user.phi.shape[0]  # number of integration points
    ndf = user.phi.shape[1]  # number of velocity dofs per spatial direction
    ndfp = user.psi.shape[1]  # number of degrees of freedom of the pressure

    # compute mapping of reference to real element
    _F, Finv, detF = isoparametric_deformation(coor, user.dphi)

    # position of the integration points
    xg = user.phi @ coor

    if user.coorsys == 1:
        # axisymmetric
        detF = 2 * np.pi * xg[:, 1] * detF

    # compute derivative of the basis functions with respect to the real
    # coordinates
    dphidx = np.zeros((ninti, ndf, ndim))

    for ip in range(ninti):
        dphidx[ip, :, :] = user.dphi[ip, :, :] @ Finv[ip, :, :]

    # pointers in unknowns
    i1 = ndf
    i2 = 2*ndf
    i3 = 2*ndf + ndfp

    # compute element matrix
    Suu = np.zeros((ndf, ndf))
    Svv = np.zeros((ndf, ndf))
    Suv = np.zeros((ndf, ndf))
    Lu = np.zeros((ndfp, ndf))
    Lv = np.zeros((ndfp, ndf))
    work = np.zeros(ninti)

    # diagonal blocks
    for N in range(ndf):
        for M in range(N, ndf):
            work[:] = 0
            for ip in range(ninti):
                for k in range(ndim):
                    work[ip] += dphidx[ip, N, k] * dphidx[ip, M, k]
            work1 = work + dphidx[:, N, 0] * dphidx[:, M, 0]
            Suu[N, M] = user.mu * np.sum(work1 * detF * user.wg)
            Suu[M, N] = Suu[N, M]  # symmetry
            work1 = work + dphidx[:, N, 1] * dphidx[:, M, 1]
            if user.coorsys == 1:
                work1 += 2 * user.phi[:, N] * user.phi[:, M] / xg[:, 1]**2
            Svv[N, M] = user.mu * np.sum(work1 * detF * user.wg)
            Svv[M, N] = Svv[N, M]  # symmetry

    # off-diagonal blocks
    for N in range(ndf):
        for M in range(ndf):
            work = dphidx[:, N, 1] * dphidx[:, M, 0]
            Suv[N, M] = user.mu * np.sum(work * detF * user.wg)

    # velocity pressure blocks
    for N in range(ndfp):
        for M in range(ndf):
            Lu[N, M] = np.sum(user.psi[:, N] * dphidx[:, M, 0] * detF
                              * user.wg)
            if user.coorsys == 1:
                work = dphidx[:, M, 1] + user.phi[:, M] / xg[:, 1]
                Lv[N, M] = np.sum(user.psi[:, N] * work * detF * user.wg)
            else:
                Lv[N, M] = np.sum(user.psi[:, N] * dphidx[:, M, 1] * detF
                                  * user.wg)

    elemmat = np.zeros((2*ndf + ndfp, 2*ndf + ndfp))
    elemmat[0:i1, 0:i1] = Suu
    elemmat[0:i1, i1:i2] = Suv
    elemmat[i1:i2, 0:i1] = Suv.T
    elemmat[i1:i2, i1:i2] = Svv
    elemmat[i2:i3, 0:i1] = -Lu
    elemmat[i2:i3, i1:i2] = -Lv
    elemmat[0:i1, i2:i3] = -Lu.T
    elemmat[i1:i2, i2:i3] = -Lv.T
    elemmat[i2:i3, i2:i3] = 0

    # compute element vector
    elemvec = np.zeros(2*ndf + ndfp)

    if user.funcnr > 0:
        fg = np.zeros((ninti, ndim))
        for ip in range(ninti):
            fg[ip, :] = user.func(user.funcnr, xg[ip, :])
        tmp = np.zeros((ninti, ndim))
        for j in range(ndim):
            for N in range(ndf):
                tmp[N, j] = np.sum(fg[:, j] * user.phi[:, N] * detF * user.wg)
        elemvec[0:i2] = tmp.reshape(ndim*ndf)

    return elemmat, elemvec


def stokes_deriv(elem: int, coor: np.ndarray, user: "User", pos: list[list[int]]):
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
    _F, Finv, _detF = isoparametric_deformation(coor, user.dphi)

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


def stokes_natboun_curve(elem: int, coor: np.ndarray, user: "User", pos: list[list[int]]):
    """
    Compute the boundary element for a natural boundary on a curve for the
    Stokes equation.

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
        The element matrix.
    elemvec : ndarray
        The element vector.

    Notes
    -----
    This function must be called in `add_boundary_elements` using
    `posvectors=0` (default).

    """

    # set some values
    ninti = user.phi.shape[0]  # number of integration points
    ndf = user.phi.shape[1]  # Number of velocity dofs per spatial direction
    ndim = coor.shape[1]  # dimension of space

    # compute mapping of reference to real element
    _dxdxi, curvel, _normal = isoparametric_deformation_curve(coor, user.dphi)

    # position of the integration points
    xg = np.dot(user.phi, coor)

    if user.coorsys == 1:
        # axisymmetric
        curvel = 2 * np.pi * xg[:, 1] * curvel

    # compute element vector
    elemvec = np.zeros(ndim * ndf)

    if user.funcnr > 0:
        fg = np.zeros((ninti, ndim))

        for ip in range(ninti):
            fg[ip, :] = user.func(user.funcnr, xg[ip, :])

        tmp = np.zeros((ndf, ndim))

        for j in range(ndim):
            for N in range(ndf):
                tmp[N, j] = np.sum(fg[:, j] * user.phi[:, N] * curvel
                                   * user.wg)

        elemvec = tmp.reshape(ndim * ndf, order='F')

    return elemvec


def stokes_flowrate_curve(elem: int, coor: np.ndarray, user: "User", pos: list[list[int]]):
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

    # Set some values
    ninti = user.phi.shape[0]  # Number of integration points
    ndf = user.phi.shape[1]    # Number of velocity dofs per spatial direction
    ndim = coor.shape[1]  # Dimension of space

    # Compute mapping of reference to real element
    _dxdxi, curvel, normal = isoparametric_deformation_curve(coor, user.dphi)

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


def stokes_pressure(elem: int, coor: np.ndarray, user: "User", pos: list[list[int]]):
    """
    Compute the pressure for post-processing in Stokes flow.

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
        The pressure at all nodes.

    Notes
    -----
    This function is intended for use with `deriv_vector` or `plot_sol`.

    """

    elemvec = user.psi @ user.u[pos[1]]

    return elemvec
