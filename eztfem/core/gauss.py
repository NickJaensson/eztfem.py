'''
Module to compute Gauss-Legendre integration points and weights.
'''
from pathlib import Path
import numpy as np


def gauss_legendre(shape, **kwargs):
    """
    Determine Gauss-Legendre integration points and weights in elements.

    Parameters
    ----------
    shape : {'line', 'quad', 'triangle'}
        Shape of the domain:
        - 'line': interval [-1, 1]
        - 'quad': domain [-1, 1] x [-1, 1]
        - 'triangle': reference triangle, left-lower half of [0, 1] x [0, 1]

    Keyword arguments
    -----------------
    n : int, optional
        The number of integration points in one direction. This only
        applies to shape='line' and shape='quad', where the number of
        integration points will be n and n^2, respectively.
    p : int, optional
        The order of the integration rule (order polynomial integrated
        exact). This only applies to shape='triangle'.

    Returns
    -------
    x : numpy.ndarray of shape (ninti,ndim)
        Coordinates of the integration points, where ndim is the spatial
        dimension and ninti is the number of integration points.
    w : numpy.ndarray of shape (ninti,)
        Weights of the integration scheme.

    """

    n = kwargs.get('n', -1)
    p = kwargs.get('p', -1)

    if shape == 'line':
        if n < 0:
            raise ValueError('n must be specified for the integration rule')
        x, w = _gauss_legendre_line(n)
    elif shape == 'quad':
        if n < 0:
            raise ValueError('n must be specified for the integration rule')
        x, w = _gauss_legendre_quad(n)
    elif shape == 'triangle':
        if p < 0:
            raise ValueError('p must be specified for the integration rule')
        x, w = _gauss_legendre_triangle(p)
    else:
        if isinstance(shape, str):
            ch = f'Invalid shape: {shape}'
        else:
            ch = 'shape must be a string'
        raise ValueError(ch)

    return x, w


def _gauss_legendre_line(n):
    """
    Compute Gauss-Legendre integration points and weights for a line segment.
    Gauss Legendre in 1D defined on the interval [-1,1]

    Parameters
    ----------
    n : int
        Number of integration points.

    Returns
    -------
    x : numpy.ndarray
        Coordinates of integration points in [-1, 1].
    w : numpy.ndarray
        Weights of the integration points.

    """
    if not isinstance(n, (int, np.integer)):
        raise TypeError(f"n must be an integer, got {type(n).__name__}")
    if n < 1:
        raise ValueError(f"n must be >= 1, got {n}")

    x, w = np.polynomial.legendre.leggauss(n)

    return x, w


def _gauss_legendre_quad(n):
    """
    Compute Gauss-Legendre integration points and weights for a quadrilateral.
    Gauss Legendre in 2D defined on the region [-1,1]x[-1,1]

    Parameters
    ----------
    n : int
        Number of integration points in one dimension.

    Returns
    -------
    x : numpy.ndarray
        Reference coordinates of the integration points, shape (n^2, 2).
    w : numpy.ndarray
        Weights of the integration points, shape (n^2,).

    """

    x1, w1 = _gauss_legendre_line(n)

    m = n**2
    x = np.zeros((m, 2))
    w = np.zeros(m)

    for i in range(n):
        for j in range(n):
            ip = i + n * (j)
            x[ip, 0] = x1[i]
            x[ip, 1] = x1[j]
            w[ip] = w1[i] * w1[j]

    return x, w


def _gauss_legendre_triangle(p):
    """
    Compute Gauss-Legendre integration points and weights for a triangle.
    Gauss Legendre in 2D defined triangle, left-lower half of [0, 1] x [0, 1]

    Parameters
    ----------
    p : int
        Order of the integration rule (order polynomial integrated exact).

    Returns
    -------
    x : numpy.ndarray
        Reference coordinates of the integration points, shape (ni, 2).
    w : numpy.ndarray
        Weights of the integration points, shape (ni,).

    Notes
    -----
    The reference region of a triangle is defined as:

    ::

            1 |\\               eta ∈ [0, 1]
            ^ |  \\             xi  ∈ [0, 1]
        eta | |    \\           xi + eta <= 1
              |      \\
            0 ---------
              0  xi -> 1

    The numerical values are based on:

    L. Zhang, T. Cui, H. Liu, "A set of symmetric quadrature rules on triangles
    and tetrahedra", Journal of Computational Mathematics, Vol. 27, No. 1,
    2009, 89-96.

    The numerical values are stored in CSV files in the
    'gauss_legendre_triangle' folder.

    """

    # number of integration points
    ni = np.array([1,  3,  6,  6,  7, 12, 15, 16, 19, 25,
                  28, 33, 37, 46, 52, 55, 61, 72, 73, 88, 91])

    x = np.zeros((ni[p-1], 2))
    w = np.zeros(ni[p-1])

    if p < 1 or p > 21:
        raise ValueError('p must be in the range [1, 21]')

    data_dir = Path(__file__).resolve().parent / "gauss_legendre_triangle"
    csv_path = data_dir / f"gauss_legendre_triangle_p{p:02d}.csv"
    data = np.loadtxt(csv_path, delimiter=",", ndmin=2, skiprows=1,
                      dtype=float)

    if not csv_path.exists():
        raise FileNotFoundError(f"Missing quadrature file: {csv_path}")

    x = data[:, :2]
    w = data[:, 2]

    return x, w
