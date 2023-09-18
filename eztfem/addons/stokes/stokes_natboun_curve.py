import numpy as np
from ...src.isoparametric_deformation_curve import isoparametric_deformation_curve

def stokes_natboun_curve(elem, coor, user, pos):
    # set some values
    ninti = user.phi.shape[0]  # number of integration points
    ndf = user.phi.shape[1]  # number of degrees of freedom of the velocity for each spatial direction
    ndim = coor.shape[1]  # dimension of space

    # compute mapping of reference to real element
    dxdxi, curvel, normal = isoparametric_deformation_curve(coor, user.dphi)

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
                tmp[N, j] = np.sum(fg[:, j] * user.phi[:, N] * curvel * user.wg)

        elemvec = tmp.reshape(ndim * ndf, order='F')

    return elemvec