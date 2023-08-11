import numpy as np

def basis_function(shape, intpol, xr):
    """
    basis function and the derivative with respect to the reference coordinates.

    Parameters:
    - shape: Shape of the element. 
             Possible values: 'line', 'quad', 'triangle'.
    - intpol: Interpolation on the element.
              Possible values for P-family: 'P0', 'P1', 'P1+', 'P2', 'P2+'
              Possible values for Q-family: 'Q1', 'Q1+', 'Q2'
              (Note: Not all combinations of shape-intpol are possible!)
    - xr: Reference coordinates where phi and dphi must be computed.
          xr[i][j] with i the point in space and j the direction in space.

    Returns:
    - phi: basis function phi[i][j], with i the point in space and j the unknown.
    - dphi: derivative of the basis function dphi[i][j][k], with i the point in space, 
            j the unknown and k the direction in space.
    """

    if shape == 'line':
        if intpol == 'P0':
            phi, dphi = basis_line_P0(xr)
        elif intpol == 'P1':
            phi, dphi = basis_line_P1(xr)
        elif intpol == 'P2':
            phi, dphi = basis_line_P2(xr)
        else:
            raise ValueError(f"Invalid intpol for shape = 'line' : {intpol}")

    elif shape == 'quad':
        if intpol == 'P0':
            phi, dphi = basis_quad_P0(xr)
        elif intpol == 'P1':
            phi, dphi = basis_quad_P1(xr)
        elif intpol == 'Q1':
            phi, dphi = basis_quad_Q1(xr)
        elif intpol == 'Q1+':
            phi, dphi = basis_quad_Q1plus(xr)
        elif intpol == 'Q2':
            phi, dphi = basis_quad_Q2(xr)
        else:
            raise ValueError(f"Invalid intpol for shape = 'quad': {intpol}")

    elif shape == 'triangle':
        if intpol == 'P0':
            phi, dphi = basis_triangle_P0(xr)
        elif intpol == 'P1':
            phi, dphi = basis_triangle_P1(xr)
        elif intpol == 'P1+':
            phi, dphi = basis_triangle_P1plus(xr)
        elif intpol == 'P2':
            phi, dphi = basis_triangle_P2(xr)
        elif intpol == 'P2+':
            phi, dphi = basis_triangle_P2plus(xr)
        else:
            raise ValueError(f"Invalid intpol for shape = 'triangle': {intpol}")

    else:
        raise ValueError(f"Invalid shape: {shape}")

    return phi, dphi

def basis_line_P0(xr):
    ni = xr.shape[0]
    nn = 1
    phi = np.zeros((ni, nn))
    dphi = np.zeros((ni, nn, 1))
    
    phi[:, :] = 1
    dphi[:, :, :] = 0
    
    return phi, dphi

def basis_line_P1(xr):
    ni = xr.shape[0]
    nn = 2
    phi = np.zeros((ni, nn))
    dphi = np.zeros((ni, nn, 1))
    
    phi[:, 0] = (1 - xr) / 2
    phi[:, 1] = (1 + xr) / 2

    dphi[:, 0, 0] = -0.5
    dphi[:, 1, 0] =  0.5

    return phi, dphi

import numpy as np

def basis_line_P2(xr):
    ni = xr.shape[0]
    nn = 3
    phi = np.zeros((ni, nn))
    dphi = np.zeros((ni, nn, 1))
    
    phi[:, 0] = -(1 - xr[:, 0]) * xr[:, 0] / 2
    phi[:, 1] = 1 - xr[:, 0] ** 2
    phi[:, 2] = (1 + xr[:, 0]) * xr[:, 0] / 2

    dphi[:, 0, 0] = -(1 - 2 * xr[:, 0]) / 2
    dphi[:, 1, 0] = -2 * xr[:, 0]
    dphi[:, 2, 0] = (1 + 2 * xr[:, 0]) / 2

    return phi, dphi

def basis_quad_P0(xr):
    ni = xr.shape[0]
    nn = 1
    phi = np.zeros((ni, nn))
    dphi = np.zeros((ni, nn, 2))
    
    phi[:, :] = 1
    dphi[:, :, :] = 0

    return phi, dphi

def basis_quad_P1(xr):
    ni = xr.shape[0]
    nn = 3
    phi = np.zeros((ni, nn))
    dphi = np.zeros((ni, nn, 2))
    
    phi[:, 0] = 1
    phi[:, 1] = xr[:, 0]
    phi[:, 2] = xr[:, 1]
    
    dphi[:, 0, :] = 0
    dphi[:, 1, 0] = 1
    dphi[:, 1, 1] = 0
    dphi[:, 2, 0] = 0
    dphi[:, 2, 1] = 1

    return phi, dphi

def basis_quad_Q1(xr):
    ni = xr.shape[0]
    nn = 4
    phi = np.zeros((ni, nn))
    dphi = np.zeros((ni, nn, 2))
    p = np.array([[1, 2], [4, 3]]).T

    phi1, dphi1 = basis_line_P1(xr[:, 0])
    phi2, dphi2 = basis_line_P1(xr[:, 1])

    for i in range(2):
        for j in range(2):
            phi[:, p[i, j]-1] = phi1[:, i] * phi2[:, j]
            dphi[:, p[i, j]-1, 0] = dphi1[:, i] * phi2[:, j]
            dphi[:, p[i, j]-1, 2] = phi1[:, i] * dphi2[:, j]

    return phi, dphi

def basis_quad_Q1plus(xr):
    ni = xr.shape[0]
    nn = 5
    phi = np.zeros((ni, nn))
    dphi = np.zeros((ni, nn, 2))
    dbubble = np.zeros((ni, 1, 2))

    phiQ1, dphiQ1 = basis_quad_Q1(xr)

    bubble = phiQ1[:, 0] * phiQ1[:, 2]

    phi[:, 0:4] = phiQ1[:, 0:4] - 4 * bubble[:, np.newaxis]
    phi[:, 4] = 16 * bubble

    for dm in range(2):
        dbubble[:, 0, dm] = (dphiQ1[:, 0, dm] * phiQ1[:, 2] +
                             phiQ1[:, 0] * dphiQ1[:, 2, dm])

    dphi[:, 0:4, :] = dphiQ1[:, 0:4, :] - 4 * dbubble
    dphi[:, 4, :] = 16 * dbubble

    return phi, dphi

def basis_quad_Q2(xr):
    ni = xr.shape[0]
    nn = 9
    phi = np.zeros((ni, nn))
    dphi = np.zeros((ni, nn, 2))
    p = np.array([[1, 2, 3], [8, 9, 4], [7, 6, 5]]).T

    phi1, dphi1 = basis_line_P2(xr[:, 0])
    phi2, dphi2 = basis_line_P2(xr[:, 1])

    for i in range(3):
        for j in range(3):
            phi[:, p[i, j]-1] = phi1[:, i] * phi2[:, j]
            dphi[:, p[i, j]-1, 0] = dphi1[:, i] * phi2[:, j]
            dphi[:, p[i, j]-1, 1] = phi1[:, i] * dphi2[:, j]

    return phi, dphi