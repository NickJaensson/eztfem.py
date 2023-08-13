import numpy as np
from src.isoparametric_deformation import isoparametric_deformation

def stokes_elem(elem, coor, user, pos):
    # STOKES_ELEM  Element routines for the Stokes equation

    # - nabla.( mu (nabla u+nabla u^T) ) + nabla p = f
    #   nabla.u = 0

    ndim = coor.shape[1]  # dimension of space
    ninti = user.phi.shape[0]  # number of integration points
    ndf = user.phi.shape[1]  # number of degrees of freedom of the velocity for each spatial direction
    ndfp = user.psi.shape[1]  # number of degrees of freedom of the pressure

    # compute mapping of reference to real element
    F, Finv, detF = isoparametric_deformation(coor, user.dphi)

    # position of the integration points
    xg = user.phi @ coor

    if user.coorsys == 1:
        # axisymmetric
        detF = 2 * np.pi * xg[:, 1] * detF

    # compute derivative of the basis functions with respect to the real coordinates
    dphidx = np.zeros((ninti, ndf, ndim))
    ndr = Finv.shape[1]

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
            Lu[N, M] = np.sum(user.psi[:, N] * dphidx[:, M, 0] * detF * user.wg)
            if user.coorsys == 1:
                work = dphidx[:, M, 1] + user.phi[:, M] / xg[:, 1]
                Lv[N, M] = np.sum(user.psi[:, N] * work * detF * user.wg)
            else:
                Lv[N, M] = np.sum(user.psi[:, N] * dphidx[:, M, 1] * detF * user.wg)

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