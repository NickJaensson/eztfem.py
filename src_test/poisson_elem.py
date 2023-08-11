import numpy as np
from src_test.isoparametric_deformation import isoparametric_deformation

def poisson_elem(elem, coor, user, pos):
    """
    Element routine for the Poisson equation: - alpha nabla^2 u = f
    
    Parameters:
    - elem: element number
    - coor: coordinates of the nodes of the element
    - user: object containing user data
    - pos: list of the positions of the degrees of freedom of each physq

    Returns:
    - elemmat: element matrix
    - elemvec: element vector
    """
    
    ninti = user.phi.shape[0]  # number of integration points
    ndf = user.phi.shape[1]    # number of degrees of freedom
    ndim = coor.shape[1]          # dimension of space

    # Compute mapping of reference to real element
    F, Finv, detF = isoparametric_deformation(coor, user.dphi)

    # Position of the integration points
    xg = np.dot(user.phi, coor)

    if user.coorsys == 1:
        # axisymmetric
        detF = 2 * np.pi * xg[:, 1] * detF

    # Compute derivative of the basis functions wrt real coordinates
    dphidx = np.zeros((ninti, ndf, ndim))
    ndr = Finv.shape[1]

    for ip in range(ninti):
        dphidx[ip] = np.reshape(user.dphi[ip], (ndf, ndr)) @ np.reshape(Finv[ip], (ndr, ndim))

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