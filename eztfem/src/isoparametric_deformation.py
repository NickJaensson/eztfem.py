import numpy as np

def isoparametric_deformation(x, dphi):
    """
    ISOPARAMETRIC_DEFORMATION  Isoparametric deformation of an element.
      [ F, Finv, detF ] = ISOPARAMETRIC_DEFORMATION ( x, dphi )
      input:
        x: coordinates of the nodes = position of the unknowns for phi and dphi
           x(i,j) with i the point in space and j the direction in space
        dphi: derivatives of the shape function with respect to the reference
              coordinates dphi(i,j,k), with i the point in space j the unknown
              and k the direction in space.
      output:
        F: Deformation gradient matrix between the reference element and the
           actual element. F(i,j,m) means that in point i the transformation is
    
                             F(j,m) = d x  / d xi
                                         j       m
        Finv: The inverse of F 
        detF: The determinant of F (Jacobian) detF(i) gives the determinant
              in point i
    """
    
    npts, ndim = dphi.shape[0], x.shape[1]
    
    F = np.zeros((npts, ndim, ndim))
    Finv = np.zeros((npts, ndim, ndim))
    detF = np.zeros(npts)

    # compute F
    for j in range(ndim):
        F[:, :, j] = dphi[:, :, j] @ x

    # compute detF and Finv
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
