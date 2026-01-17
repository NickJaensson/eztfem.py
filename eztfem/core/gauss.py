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
    num_int_points : int, optional
        The number of integration points in one direction. This only
        applies to shape='line' and shape='quad', where the number of
        integration points will be num_int_points and num_int_points^2,
        respectively.
    integration_order : int, optional
        The order of the integration rule (order polynomial integrated
        exact). This only applies to shape='triangle'.
    n : int, optional
        Deprecated. Use num_int_points instead.
    p : int, optional
        Deprecated. Use integration_order instead.

    Returns
    -------
    x : numpy.ndarray of shape (ninti,ndim)
        Coordinates of the integration points, where ndim is the spatial
        dimension and ninti is the number of integration points.
    w : numpy.ndarray of shape (ninti,)
        Weights of the integration scheme.

    """
    shape_handlers = {
        'line': ('num_int_points', _gauss_legendre_line),
        'quad': ('num_int_points', _gauss_legendre_quad),
        'triangle': ('integration_order', _gauss_legendre_triangle),
    }

    if shape not in shape_handlers:
        valid_shapes = list(shape_handlers.keys())
        raise ValueError(
            f"Invalid shape: '{shape}'. Must be one of {valid_shapes}"
        )

    param_name, handler = shape_handlers[shape]

    # Support both new and old parameter names for backward compatibility
    new_param_name = param_name
    old_param_map = {'num_int_points': 'n', 'integration_order': 'p'}
    old_param_name = old_param_map.get(param_name)

    param_value = kwargs.get(new_param_name, -1)
    if param_value < 0 and old_param_name:
        param_value = kwargs.get(old_param_name, -1)

    if param_value < 0:
        msg = f"{new_param_name} must be specified for the integration rule"
        raise ValueError(msg)

    return handler(param_value)


def _gauss_legendre_line(num_int_points):
    """
    Compute Gauss-Legendre integration points and weights for a line segment.
    Gauss Legendre in 1D defined on the interval [-1,1]

    Parameters
    ----------
    num_int_points : int
        Number of integration points.

    Returns
    -------
    coords : numpy.ndarray
        Coordinates of integration points in [-1, 1].
    weights : numpy.ndarray
        Weights of the integration points.

    """
    if not isinstance(num_int_points, (int, np.integer)):
        raise TypeError(
            f"num_int_points must be an integer, got "
            f"{type(num_int_points).__name__}"
        )
    if num_int_points < 1:
        raise ValueError(f"num_int_points must be >= 1, got {num_int_points}")

    coords_1d, weights = np.polynomial.legendre.leggauss(num_int_points)

    return coords_1d, weights


def _gauss_legendre_quad(num_int_points):
    """
    Compute Gauss-Legendre integration points and weights for a quadrilateral.
    Gauss Legendre in 2D defined on the region [-1,1]x[-1,1]

    Parameters
    ----------
    num_int_points : int
        Number of integration points in one dimension.

    Returns
    -------
    coords_2d : numpy.ndarray
        Reference coordinates of the integration points,
        shape (num_int_points^2, 2).
    weights_2d : numpy.ndarray
        Weights of the integration points, shape (num_int_points^2,).

    """
    line_points, line_weights = _gauss_legendre_line(num_int_points)

    # Create meshgrid of integration points with C-order (row-major)
    xi, eta = np.meshgrid(line_points, line_points, indexing='xy')
    coords_2d = np.column_stack([xi.ravel(), eta.ravel()])

    # Compute tensor product of weights
    weights_grid = np.outer(line_weights, line_weights)
    weights_2d = weights_grid.ravel()

    return coords_2d, weights_2d


def _gauss_legendre_triangle(integration_order):
    """
    Compute Gauss-Legendre integration points and weights for a triangle.
    Gauss Legendre in 2D defined triangle, left-lower half of [0, 1] x [0, 1]

    Parameters
    ----------
    integration_order : int
        Order of the integration rule (order polynomial integrated exact).
        Must be in the range [1, 21].

    Returns
    -------
    coords_tri : numpy.ndarray
        Reference coordinates of the integration points, shape (ni, 2).
    weights_tri : numpy.ndarray
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
    if 1 <= integration_order <= 21:
        quadrature_dir = (
            Path(__file__).resolve().parent / "gauss_legendre_triangle"
        )
        quadrature_file = (
            quadrature_dir /
            f"gauss_legendre_triangle_p{integration_order:02d}.csv"
        )

        if not quadrature_file.exists():
            raise FileNotFoundError(
                f"Missing quadrature file: {quadrature_file}"
            )

        quadrature_data = np.loadtxt(
            quadrature_file, delimiter=",", ndmin=2, skiprows=1, dtype=float
        )

        coords_tri = quadrature_data[:, :2]
        weights_tri = quadrature_data[:, 2]
        return coords_tri, weights_tri

    raise ValueError('order must be in the range [1, 21]')
