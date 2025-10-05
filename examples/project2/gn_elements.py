import numpy as np
from eztfem import isoparametric_deformation


def viscosity_function(gammadot, user):
    """
    Calculate the non-Newtonian viscosity for a given shear rate

    Parameters
    ----------
    gammadot : float
        Shear rate.
    user : User
        User object containing shape function, Gauss points, parameters
        user.u the solution vector from the previous iteration
        user.gnmodel the generalized Newtonian model:
              1: power-law, 2: Carreau, 3: Carreau-Yasuda
        user.m factor m of the Power-law model (for gnmodel=1)
        user.n Exponent n (for gnmodel=1, 2 and 3)
        user.eta0 zero-shear viscosity (for gnmodel=2, 3)
        user.etainf infinite-shear viscosity (for gnmodel=2, 3)
        user.lambda lambda parameter (for gnmodel=2, 3)
        user.a parameter a for the Carreau-Yasuda model (gnmodel=3)

    Returns
    -------
    eta : float
        The viscosity.

    """

    if user.gnmodel == 1:
        # Power-law
        eta = user.m / gammadot ** (1 - user.n)
    elif user.gnmodel == 2:
        # Carreau
        eta = (user.eta_inf
               + (user.eta_0 - user.eta_inf)
               / (1 + (user.llambda_ * gammadot) ** 2) ** ((1 - user.n) / 2))
    elif user.gnmodel == 3:
        # Carreau-Yasuda
        eta = (user.eta_inf +
               (user.eta_0 - user.eta_inf)
               / (1 + (user.llambda_ * gammadot) ** user.a)
               ** ((1 - user.n) / user.a))
    else:
        raise ValueError("Invalid gnmodel value")

    return eta


def gn_elem(elem, coor, user, pos):
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
        user.u the solution vector from the previous iteration
        user.gnmodel the generalized Newtonian model:
              1: power-law, 2: Carreau, 3: Carreau-Yasuda
        user.m factor m of the Power-law model (for gnmodel=1)
        user.n Exponent n (for gnmodel=1, 2 and 3)
        user.eta0 zero-shear viscosity (for gnmodel=2, 3)
        user.etainf infinite-shear viscosity (for gnmodel=2, 3)
        user.lambda lambda parameter (for gnmodel=2, 3)
        user.a parameter a for the Carreau-Yasuda model (gnmodel=3)
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
    F, Finv, detF = isoparametric_deformation(coor, user.dphi)

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

    # get velocity vector from previous iteration

    u = user.u[pos[0]]  # velocity

    uvector = np.reshape(u, (ndf, ndim), order='F')

    # compute velocity gradient tensor

    gradu = np.zeros((ninti, ndim, ndim))

    for j in range(ndim):
        gradu[:, :, j] = dphidx[:, :, j] @ uvector

    # compute effective strain rate in all integration points

    gd = np.zeros(ninti)
    for ip in range(ninti):
        # second invariant
        gu = np.reshape(gradu[ip, :, :], (ndim, ndim))
        D = (gu + gu.T) / 2
        gd[ip] = np.tensordot(D, D)

    gd = np.sqrt(2*gd)

    # compute viscosity eta in all integration points

    # TODO Compute the ninti-by-1 array (vector) of the viscosity function
    # Use the material parameters in the user structure

    eta = viscosity_function(gd, user)

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
            Suu[N, M] = np.sum(eta * work1 * detF * user.wg)
            Suu[M, N] = Suu[N, M]  # symmetry
            work1 = work + dphidx[:, N, 1] * dphidx[:, M, 1]
            if user.coorsys == 1:
                work1 += 2 * user.phi[:, N] * user.phi[:, M] / xg[:, 1]**2
            Svv[N, M] = np.sum(eta * work1 * detF * user.wg)
            Svv[M, N] = Svv[N, M]  # symmetry

    # off-diagonal blocks
    for N in range(ndf):
        for M in range(ndf):
            work = dphidx[:, N, 1] * dphidx[:, M, 0]
            Suv[N, M] = np.sum(eta * work * detF * user.wg)

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
