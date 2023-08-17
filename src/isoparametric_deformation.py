import numpy as np

def isoparametric_deformation(x, dphi):
    """
    Isoparametric deformation of an element.
    
    Parameters:
    - x: coordinates of the nodes
    - dphi: derivatives of the shape function wrt reference coordinates
    
    Returns:
    - F: Deformation gradient matrix
    - Finv: Inverse of F
    - detF: Determinant of F
    """
    
    npts, ndim = dphi.shape[0], x.shape[1]
    
    F = np.zeros((npts, ndim, ndim))
    Finv = np.zeros((npts, ndim, ndim))
    detF = np.zeros(npts)

    # Compute F
    for j in range(ndim):
        F[:, :, j] = dphi[:, :, j] @ x

    # Compute detF and Finv
    if ndim == 1:
        # 1D
        for ip in range(npts):
            detF[ip] = F[ip, 0, 0]
            Finv[ip, 0, 0] = 1 / F[ip, 0, 0]
    elif ndim == 2:
        # 2D
        for ip in range(npts):
            detF[ip] = F[ip, 0, 0] * F[ip, 1, 1] - F[ip, 0, 1] * F[ip, 1, 0]
            Finv[ip, 0, 0] = F[ip, 1, 1] / detF[ip]
            Finv[ip, 0, 1] = -F[ip, 0, 1] / detF[ip]
            Finv[ip, 1, 0] = -F[ip, 1, 0] / detF[ip]
            Finv[ip, 1, 1] = F[ip, 0, 0] / detF[ip]

    return F, Finv, detF
