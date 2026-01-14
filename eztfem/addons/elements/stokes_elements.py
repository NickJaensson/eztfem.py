'''
Element routines for the Stokes equation.
'''
import numpy as np
from ...core.shapefunc import isoparametric_deformation, \
    isoparametric_deformation_curve


def stokes_elem(elem, coor, user, pos):
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
    _, fmat_inv, det_fmat = isoparametric_deformation(coor, user.dphi)

    # position of the integration points
    xg = user.phi @ coor

    if user.coorsys == 1:
        # axisymmetric
        det_fmat = 2 * np.pi * xg[:, 1] * det_fmat

    # compute derivative of the basis functions with respect to the real
    # coordinates
    dphidx = np.zeros((ninti, ndf, ndim))

    for ip in range(ninti):
        dphidx[ip, :, :] = user.dphi[ip, :, :] @ fmat_inv[ip, :, :]

    # pointers in unknowns
    i1 = ndf
    i2 = 2*ndf
    i3 = 2*ndf + ndfp

    # compute element matrix
    smat_uu = np.zeros((ndf, ndf))
    smat_vv = np.zeros((ndf, ndf))
    smat_uv = np.zeros((ndf, ndf))
    lmat_u = np.zeros((ndfp, ndf))
    lmat_v = np.zeros((ndfp, ndf))
    work = np.zeros(ninti)

    # diagonal blocks
    for idf in range(ndf):
        for jdf in range(idf, ndf):
            work[:] = 0
            for ip in range(ninti):
                for k in range(ndim):
                    work[ip] += dphidx[ip, idf, k] * dphidx[ip, jdf, k]
            work1 = work + dphidx[:, idf, 0] * dphidx[:, jdf, 0]
            smat_uu[idf, jdf] = user.mu * np.sum(work1 * det_fmat * user.wg)
            smat_uu[jdf, idf] = smat_uu[idf, jdf]  # symmetry
            work1 = work + dphidx[:, idf, 1] * dphidx[:, jdf, 1]
            if user.coorsys == 1:
                work1 += 2 * user.phi[:, idf] * user.phi[:, jdf] / xg[:, 1]**2
            smat_vv[idf, jdf] = user.mu * np.sum(work1 * det_fmat * user.wg)
            smat_vv[jdf, idf] = smat_vv[idf, jdf]  # symmetry

    # off-diagonal blocks
    for idf in range(ndf):
        for jdf in range(ndf):
            work = dphidx[:, idf, 1] * dphidx[:, jdf, 0]
            smat_uv[idf, jdf] = user.mu * np.sum(work * det_fmat * user.wg)

    # velocity pressure blocks
    for idf in range(ndfp):
        for jdf in range(ndf):
            lmat_u[idf, jdf] = np.sum(user.psi[:, idf] * dphidx[:, jdf, 0]
                                      * det_fmat * user.wg)
            if user.coorsys == 1:
                work = dphidx[:, jdf, 1] + user.phi[:, jdf] / xg[:, 1]
                lmat_v[idf, jdf] = np.sum(user.psi[:, idf] * work * det_fmat
                                          * user.wg)
            else:
                lmat_v[idf, jdf] = np.sum(user.psi[:, idf] * dphidx[:, jdf, 1]
                                          * det_fmat * user.wg)

    elemmat = np.zeros((2*ndf + ndfp, 2*ndf + ndfp))
    elemmat[0:i1, 0:i1] = smat_uu
    elemmat[0:i1, i1:i2] = smat_uv
    elemmat[i1:i2, 0:i1] = smat_uv.T
    elemmat[i1:i2, i1:i2] = smat_vv
    elemmat[i2:i3, 0:i1] = -lmat_u
    elemmat[i2:i3, i1:i2] = -lmat_v
    elemmat[0:i1, i2:i3] = -lmat_u.T
    elemmat[i1:i2, i2:i3] = -lmat_v.T
    elemmat[i2:i3, i2:i3] = 0

    # compute element vector
    elemvec = np.zeros(2*ndf + ndfp)

    if user.funcnr > 0:
        fg = np.zeros((ninti, ndim))
        for ip in range(ninti):
            fg[ip, :] = user.func(user.funcnr, xg[ip, :])
        tmp = np.zeros((ninti, ndim))
        for j in range(ndim):
            for idf in range(ndf):
                tmp[idf, j] = np.sum(fg[:, j] * user.phi[:, idf] * det_fmat
                                     * user.wg)
        elemvec[0:i2] = tmp.reshape(ndim*ndf)

    return elemmat, elemvec


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
    _, fmat_inv, _ = isoparametric_deformation(coor, user.dphi)

    # Compute derivative of the basis functions with respect to the real
    # coordinates
    dphidx = np.zeros((nodalp, ndf, ndim))

    for ip in range(nodalp):
        dphidx[ip, :, :] = user.dphi[ip, :, :].dot(fmat_inv[ip, :, :])

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


def stokes_natboun_curve(elem, coor, user, pos):
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
    _, curvel, _ = isoparametric_deformation_curve(coor, user.dphi)

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
            for idf in range(ndf):
                tmp[idf, j] = np.sum(fg[:, j] * user.phi[:, idf] * curvel
                                     * user.wg)

        elemvec = tmp.reshape(ndim * ndf, order='F')

    return elemvec


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

    # Set some values
    ninti = user.phi.shape[0]  # Number of integration points
    ndf = user.phi.shape[1]    # Number of velocity dofs per spatial direction
    ndim = coor.shape[1]  # Dimension of space

    # Compute mapping of reference to real element
    _, curvel, normal = isoparametric_deformation_curve(coor, user.dphi)

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


def stokes_pressure(elem, coor, user, pos):
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
