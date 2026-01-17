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
    shape_handlers = {
        'line': ('n', _gauss_legendre_line),
        'quad': ('n', _gauss_legendre_quad),
        'triangle': ('p', _gauss_legendre_triangle),
    }

    if shape not in shape_handlers:
        valid_shapes = list(shape_handlers.keys())
        raise ValueError(
            f"Invalid shape: '{shape}'. Must be one of {valid_shapes}"
        )

    param_name, handler = shape_handlers[shape]
    param_value = kwargs.get(param_name, -1)

    if param_value < 0:
        msg = f"{param_name} must be specified for the integration rule"
        raise ValueError(msg)

    return handler(param_value)


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

    # Create meshgrid of integration points with C-order (row-major)
    xi, eta = np.meshgrid(x1, x1, indexing='xy')
    x = np.column_stack([xi.ravel(), eta.ravel()])

    # Compute tensor product of weights
    w_2d = np.outer(w1, w1)
    w = w_2d.ravel()

    return x, w


def _gauss_legendre_triangle(p):
    """
    Compute Gauss-Legendre integration points and weights for a triangle.
    Gauss Legendre in 2D defined triangle, left-lower half of [0, 1] x [0, 1]

    Parameters
    ----------
    p : int
        Order of the integration rule (order polynomial integrated exact).
        Must be in the range [1, 21].

    Returns
    -------
    x : numpy.ndarray
        Reference coordinates of the integration points, shape (ni, 2).
    w : numpy.ndarray
        Weights of the integration points, shape (ni,).

    Notes
    -----
    Number of integration points for each order (p=order, ni=nr of points)
    p   ni     p  ni     p  ni     p  ni     p  ni
    01  01     06 12     11 28     16 55     21 91
    02  03     07 15     12 33     17 61
    03  06     08 16     13 37     18 72
    04  06     09 19     14 46     19 73
    05  07     10 25     15 52     20 88

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
    if 1 <= p <= 21:
        data_dir = Path(__file__).resolve().parent / "gauss_legendre_triangle"
        csv_path = data_dir / f"gauss_legendre_triangle_p{p:02d}.csv"

        if not csv_path.exists():
            raise FileNotFoundError(f"Missing quadrature file: {csv_path}")

        data = np.loadtxt(
            csv_path, delimiter=",", ndmin=2, skiprows=1, dtype=float
        )

        return data[:, :2], data[:, 2]

    raise ValueError('p must be in the range [1, 21]')
