import numpy as np
from eztfem.src.isoparametric_deformation_curve import isoparametric_deformation_curve

def streamfunction_natboun_curve(elem, coor, user, pos, posvec):

    # Set some values
    ninti = user.phi.shape[0]  # Number of integration points
    ndf = user.phi.shape[1]  # Number of degrees of freedom
    ndim = coor.shape[1]  # Dimension of space

    # Compute mapping of reference to real element
    dxdxi, curvel, normal = isoparametric_deformation_curve(coor, user.dphi)

    # Position of the integration points
    xg = user.phi @ coor

    if user.coorsys == 1:
        # Axisymmetric
        curvel = 2 * np.pi * xg[:, 1] * curvel

    # Compute element vector
    u = user.v[posvec[1]]  # Velocity
    tmp = u.reshape((ndf, ndim),order='F')

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
