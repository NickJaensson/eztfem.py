import numpy as np
from ...src.isoparametric_deformation import isoparametric_deformation

def streamfunction_elem(elem, coor, user, pos, posvec):
    # Set some values
    ninti = user.phi.shape[0]  # number of integration points
    ndf = user.phi.shape[1]  # number of degrees of freedom
    ndim = coor.shape[1]  # dimension of space

    # compute mapping of reference to real element
    F, Finv, detF = isoparametric_deformation(coor, user.dphi)

    # Position of the integration points
    xg = user.phi @ coor

    if user.coorsys == 1:
        # Axisymmetric
        detF = 2 * np.pi * xg[:, 1] * detF

    # compute derivative of the basis functions with respect to the real coordinates
    dphidx = np.zeros((ninti, ndf, ndim))
    ndr = Finv.shape[1]

    for ip in range(ninti):
        dphidx[ip, :, :] = user.dphi[ip, :, :] @ Finv[ip, :, :]

    # Compute element matrix
    elemmat = np.zeros((ndf, ndf))
    work = np.zeros(ninti)

    for i in range(ndf):
        for j in range(i, ndf):
            work[:] = 0
            for ip in range(ninti):
                for k in range(ndim):
                    work[ip] += dphidx[ip, i, k] * dphidx[ip, j, k]
            elemmat[i, j] = np.sum(work * detF * user.wg)
            elemmat[j, i] = elemmat[i, j]  # Symmetry

    # Compute element vector
    u = user.v[posvec[1]]  # Velocity
    tmp = u.reshape((ndf, ndim),order='F')

    dudy = np.dot(dphidx[:, :, 1], tmp[:, 0])
    dvdx = np.dot(dphidx[:, :, 0], tmp[:, 1])

    if user.coorsys == 1:
        work = np.dot(user.phi, u[:ndf])
        elemvec = 2 * np.pi * np.dot(user.phi.T, (xg[:, 1] * (dvdx - dudy) - work) * detF * user.wg)
    else:
        elemvec = np.dot(user.phi.T, (dvdx - dudy) * detF * user.wg)

    return elemmat, elemvec
