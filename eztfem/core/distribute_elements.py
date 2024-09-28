import numpy as np


def distribute_elements(nelem, ratio, factor):
    """
    Generate non-equidistant n elements on the interval [0, 1].

    Parameters
    ----------
    n : int
        Number of elements.
    ratio : int
        Determines the distribution of elements:
        0: Equidistant mesh.
        1: The size of the last element is `factor` times the first.
        2: The size of an element is `factor` times the previous one.
        3: The size of the last element is `1/factor` times the first.
        4: The size of an element is `1/factor` times the previous one.
    factor : float
        Factor used in the distribution calculations.

    Returns
    -------
    x : numpy.ndarray
        Coordinates of n+1 points.

    Notes
    -----
    The interval [0, 1] is divided into n elements:

    .. math::
        dx_{i+1} = g \\cdot dx_{i-1}

    .. math::
        1 = (1 + g + g^2 + g^3 + \\ldots + g^{n-1}) \\cdot dx_1 =
        \\frac{1 - g^n}{1 - g} \\cdot dx_1

    With :math:`fac = 1 + g + g^2 + \\ldots + g^{n-1}`, we have:

    .. math::
        dx_1 = \\frac{1}{fac} = \\frac{1 - g}{1 - g^n}

    .. math::
        dx_n = g^{n-1} \\cdot dx_1
    """

    if factor < 0:
        raise ValueError("Negative factor")

    if nelem <= 1:
        raise ValueError("Number of elements must be at least 2")

    match ratio:
        case 0:
            g = 1
        case 1:
            g = np.exp(np.log(factor) / (nelem - 1))
        case 2:
            g = factor
        case 3:
            g = np.exp(-np.log(factor) / (nelem - 1))
        case 4:
            g = 1 / factor
        case _:
            raise ValueError(f"Invalid value for ratio: {ratio}")

    # generate mesh
    fac = 1.0
    for i in range(1, nelem):
        fac = fac + g**i

    # size of first element
    dx_1 = 1.0 / fac

    # generate all elements
    x = np.zeros(nelem+1)
    x[0] = 0
    dx = dx_1
    for i in range(1, nelem+1):
        x[i] = x[i-1] + dx
        dx = g * dx

    # test whether x(n+1) = 1
    if abs(x[-1] - 1) > 1e-10:
        raise ValueError(f"End value x(n+1) != 1: {x[-1]:.5e}")
    else:
        x[-1] = 1.0

    return x
