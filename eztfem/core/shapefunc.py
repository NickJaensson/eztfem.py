import numpy as np


def basis_function(shape, intpol, xr):
    """
    Compute the basis function and its derivative at reference coordinates.

    Parameters
    ----------
    shape : {'line', 'quad', 'triangle'}
        Shape of the element:
        - 'line': line element on the reference interval [-1, 1]
        - 'quad': quadrilateral element on the reference domain
        [-1, 1] x [-1, 1]
        - 'triangle': triangular element on the reference domain left-lower
        half of [0, 1] x [0, 1]
    intpol : {'P0', 'P1', 'P1+', 'P2', 'P2+', 'Q1', 'Q1+', 'Q2'}
        Interpolation on the element:
        - P-family: 'P0', 'P1', 'P1+', 'P2', 'P2+'
        - Q-family: 'Q1', 'Q1+', 'Q2'
        Note: not all combinations of shape and intpol are possible!
    xr : array_like
        Reference coordinates where phi and dphi must be computed.
        xr[i, j] with i the point in space and j the direction in space.

    Returns
    -------
    phi : numpy.ndarray
        Basis function phi[i, j], with i the point in space and j the unknown.
    dphi : numpy.ndarray
        Derivative of the basis function dphi[i, j, k], with i the point in
        space, j the unknown, and k the direction in space.

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
            raise ValueError(f"Invalid intpol for shape = 'triangle': \
                             {intpol}")

    else:
        raise ValueError(f"Invalid shape: {shape}")

    return phi, dphi


def basis_triangle_P0(xr):
    ni = xr.shape[0]
    nn = 1
    phi = np.zeros((ni, nn))
    dphi = np.zeros((ni, nn, 2))

    phi[:, :] = 1
    dphi[:, :, :] = 0

    return phi, dphi


def basis_triangle_P1(xr):
    phi, dphi = barycentric(xr)
    return phi, dphi


def basis_triangle_P1plus(xr):

    ni, _ = xr.shape
    nn = 4
    phi = np.zeros((ni, nn))
    dphi = np.zeros((ni, nn, 2))
    dbubble = np.zeros((ni, 1, 2))

    lam, dlam = barycentric(xr)

    bubble = lam[:, 0] * lam[:, 1] * lam[:, 2]

    phi[:, 0] = lam[:, 0] - 9 * bubble
    phi[:, 1] = lam[:, 1] - 9 * bubble
    phi[:, 2] = lam[:, 2] - 9 * bubble
    phi[:, 3] = 27 * bubble

    for dm in range(2):
        dbubble[:, 0, dm] = (dlam[:, 0, dm] * lam[:, 1] * lam[:, 2] +
                             lam[:, 0] * dlam[:, 1, dm] * lam[:, 2] +
                             lam[:, 0] * lam[:, 1] * dlam[:, 2, dm])

    dphi[:, 0, :] = dlam[:, 0, :] - 9 * dbubble[:, 0, :]
    dphi[:, 1, :] = dlam[:, 1, :] - 9 * dbubble[:, 0, :]
    dphi[:, 2, :] = dlam[:, 2, :] - 9 * dbubble[:, 0, :]
    dphi[:, 3, :] = 27 * dbubble[:, 0, :]

    return phi, dphi


def basis_triangle_P2(xr):
    ni, _ = xr.shape
    nn = 6
    phi = np.zeros((ni, nn))
    dphi = np.zeros((ni, nn, 2))

    lam, dlam = barycentric(xr)

    phi[:, 0] = lam[:, 0] * (2 * lam[:, 0] - 1)
    phi[:, 1] = 4 * lam[:, 0] * lam[:, 1]
    phi[:, 2] = lam[:, 1] * (2 * lam[:, 1] - 1)
    phi[:, 3] = 4 * lam[:, 1] * lam[:, 2]
    phi[:, 4] = lam[:, 2] * (2 * lam[:, 2] - 1)
    phi[:, 5] = 4 * lam[:, 2] * lam[:, 0]

    for dm in range(2):
        dphi[:, 0, dm] = (4 * lam[:, 0] - 1) * dlam[:, 0, dm]
        dphi[:, 1, dm] = 4 * (dlam[:, 0, dm] * lam[:, 1]
                              + lam[:, 0] * dlam[:, 1, dm])
        dphi[:, 2, dm] = (4 * lam[:, 1] - 1) * dlam[:, 1, dm]
        dphi[:, 3, dm] = 4 * (dlam[:, 1, dm] * lam[:, 2]
                              + lam[:, 1] * dlam[:, 2, dm])
        dphi[:, 4, dm] = (4 * lam[:, 2] - 1) * dlam[:, 2, dm]
        dphi[:, 5, dm] = 4 * (dlam[:, 2, dm] * lam[:, 0]
                              + lam[:, 2] * dlam[:, 0, dm])

    return phi, dphi


def basis_triangle_P2plus(xr):
    ni, _ = xr.shape
    nn = 7
    phi = np.zeros((ni, nn))
    dphi = np.zeros((ni, nn, 2))
    dbubble = np.zeros((ni, 1, 2))

    lam, dlam = barycentric(xr)

    bubble = lam[:, 0] * lam[:, 1] * lam[:, 2]

    phi[:, 0] = lam[:, 0] * (2 * lam[:, 0] - 1) + 3 * bubble
    phi[:, 1] = 4 * lam[:, 0] * lam[:, 1] - 12 * bubble
    phi[:, 2] = lam[:, 1] * (2 * lam[:, 1] - 1) + 3 * bubble
    phi[:, 3] = 4 * lam[:, 1] * lam[:, 2] - 12 * bubble
    phi[:, 4] = lam[:, 2] * (2 * lam[:, 2] - 1) + 3 * bubble
    phi[:, 5] = 4 * lam[:, 2] * lam[:, 0] - 12 * bubble
    phi[:, 6] = 27 * bubble

    for dm in range(2):
        dbubble[:, 0, dm] = (dlam[:, 0, dm] * lam[:, 1] * lam[:, 2] +
                             lam[:, 0] * dlam[:, 1, dm] * lam[:, 2] +
                             lam[:, 0] * lam[:, 1] * dlam[:, 2, dm])

        dphi[:, 0, dm] = ((4 * lam[:, 0] - 1)
                          * dlam[:, 0, dm] + 3 * dbubble[:, 0, dm])
        dphi[:, 1, dm] = (4 * (dlam[:, 0, dm]
                          * lam[:, 1] + lam[:, 0] * dlam[:, 1, dm])
                          - 12 * dbubble[:, 0, dm])
        dphi[:, 2, dm] = ((4 * lam[:, 1] - 1) * dlam[:, 1, dm]
                          + 3 * dbubble[:, 0, dm])
        dphi[:, 3, dm] = (4 * (dlam[:, 1, dm] * lam[:, 2]
                          + lam[:, 1] * dlam[:, 2, dm])
                          - 12 * dbubble[:, 0, dm])
        dphi[:, 4, dm] = ((4 * lam[:, 2] - 1) * dlam[:, 2, dm]
                          + 3 * dbubble[:, 0, dm])
        dphi[:, 5, dm] = (4 * (dlam[:, 2, dm] * lam[:, 0]
                          + lam[:, 2] * dlam[:, 0, dm])
                          - 12 * dbubble[:, 0, dm])

    dphi[:, 6, :] = 27 * dbubble[:, 0, :]

    return phi, dphi


def basis_line_P0(xr):
    ni = xr.shape[0]
    nn = 1
    phi = np.zeros((ni, nn))
    dphi = np.zeros((ni, nn, 1))

    phi[:, :] = 1
    dphi[:, :] = 0

    return phi, dphi


def basis_line_P1(xr):
    ni = xr.shape[0]
    nn = 2
    phi = np.zeros((ni, nn))
    dphi = np.zeros((ni, nn, 1))

    phi[:, 0] = (1 - xr) / 2
    phi[:, 1] = (1 + xr) / 2

    dphi[:, 0, 0] = -0.5
    dphi[:, 1, 0] = 0.5

    return phi, dphi


def basis_line_P2(xr):
    ni = xr.shape[0]
    nn = 3
    phi = np.zeros((ni, nn))
    dphi = np.zeros((ni, nn, 1))

    phi[:, 0] = -(1 - xr) * xr / 2
    phi[:, 1] = 1 - xr ** 2
    phi[:, 2] = (1 + xr) * xr / 2

    dphi[:, 0, 0] = -(1 - 2 * xr) / 2
    dphi[:, 1, 0] = -2 * xr
    dphi[:, 2, 0] = (1 + 2 * xr) / 2

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
    p = np.array([[0, 1], [3, 2]]).T

    phi1, dphi1 = basis_line_P1(xr[:, 0])
    phi2, dphi2 = basis_line_P1(xr[:, 1])

    for i in range(2):
        for j in range(2):
            phi[:, p[i, j]] = phi1[:, i] * phi2[:, j]
            dphi[:, p[i, j], 0] = dphi1[:, i, 0] * phi2[:, j]
            dphi[:, p[i, j], 1] = phi1[:, i] * dphi2[:, j, 0]

    return phi, dphi


def basis_quad_Q1plus(xr):
    ni = xr.shape[0]
    nn = 5
    phi = np.zeros((ni, nn))
    dphi = np.zeros((ni, nn, 2))
    dbubble = np.zeros((ni, 2))

    phiQ1, dphiQ1 = basis_quad_Q1(xr)

    bubble = phiQ1[:, 0] * phiQ1[:, 2]

    phi[:, 0:4] = phiQ1[:, 0:4] - 4 * bubble[:, np.newaxis]
    phi[:, 4] = 16 * bubble

    for dm in range(2):
        dbubble[:, dm] = (dphiQ1[:, 0, dm] * phiQ1[:, 2] +
                          phiQ1[:, 0] * dphiQ1[:, 2, dm])

    dphi[:, 0, :] = dphiQ1[:, 0, :] - 4 * dbubble
    dphi[:, 1, :] = dphiQ1[:, 1, :] - 4 * dbubble
    dphi[:, 2, :] = dphiQ1[:, 2, :] - 4 * dbubble
    dphi[:, 3, :] = dphiQ1[:, 3, :] - 4 * dbubble

    dphi[:, 4, :] = 16 * dbubble

    return phi, dphi


def basis_quad_Q2(xr):
    ni = xr.shape[0]
    nn = 9
    phi = np.zeros((ni, nn))
    dphi = np.zeros((ni, nn, 2))
    p = np.array([[0, 1, 2], [7, 8, 3], [6, 5, 4]]).T

    phi1, dphi1 = basis_line_P2(xr[:, 0])
    phi2, dphi2 = basis_line_P2(xr[:, 1])

    for i in range(3):
        for j in range(3):
            phi[:, p[i, j]] = phi1[:, i] * phi2[:, j]
            dphi[:, p[i, j], 0] = dphi1[:, i, 0] * phi2[:, j]
            dphi[:, p[i, j], 1] = phi1[:, i] * dphi2[:, j, 0]

    return phi, dphi


def barycentric(xr):
    """
    Compute the barycentric coordinates and their derivatives for a triangular
    element at reference coordinates.

    Parameters
    ----------
    xr : array_like
        Reference coordinates where lambda and dlambda must be computed.
        xr[i, j] with i the point in space and j the direction in space.

    Returns
    -------
    lambda : numpy.ndarray
        Barycentric coordinates lambda[i, j], with i the point in space and j
        the number.
    dlambda : numpy.ndarray
        Derivative of the barycentric coordinates dlambda[i, j, k], with i the
        point in space, j the number, and k the direction in space.

    Notes
    -----

    ::

        Triangles:
        The reference region of a triangle is defined as:
                1 |\\               eta ∈ [0,1]
                ^ |  \\             xi  ∈ [0,1]
            eta | |    \\           xi + eta <= 1
                  |      \\
                0 ---------
                  0  xi -> 1
        The reference coordinates xr = (xi, eta) are in the lower triangle of the
        region [0,1] x [0,1] (xi + eta <= 1).
        The three barycentric coordinates are:
        - lambda1 = 1 - xi - eta
        - lambda2 = xi
        - lambda3 = eta

    """

    # Get sizes from xr shape
    ni, refdim = xr.shape

    # Initialize the barycentric coordinates and their derivatives
    nn = 3
    lambda_ = np.zeros((ni, nn))
    dlambda = np.zeros((ni, nn, refdim))

    # Calculate the barycentric coordinates
    lambda_[:, 0] = 1 - xr[:, 0] - xr[:, 1]
    lambda_[:, 1] = xr[:, 0]
    lambda_[:, 2] = xr[:, 1]

    # Derivatives of the barycentric coordinates
    dlambda[:, 0, 0] = -1
    dlambda[:, 0, 1] = -1
    dlambda[:, 1, 0] = 1
    dlambda[:, 1, 1] = 0
    dlambda[:, 2, 0] = 0
    dlambda[:, 2, 1] = 1

    return lambda_, dlambda
