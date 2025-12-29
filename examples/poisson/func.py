'''
Supporting functions for the Poisson example problems.
'''
import math


def func(nr, x):
    """
    Evaluate a scalar function as a function of spatial coordinates.

    Parameters
    ----------
    nr : int
        Function selector.
    x : array_like of shape (ndim,)
        Spatial coordinates, where ndim is the spatial dimension.

    Returns
    -------
    value : float
        Function value at x.
    """

    if nr == 1:
        value = (x[1] + 1) ** 2
    elif nr == 3:
        value = 1 + math.cos(math.pi * x[0]) * math.cos(math.pi * x[1])
    elif nr == 4:
        value = 2 * math.pi ** 2 * math.cos(math.pi * x[0]) \
            * math.cos(math.pi * x[1])
    elif nr == 5:
        value = 1.0
    elif nr == 6:
        value = math.cos(math.pi * x[0]) * \
            math.cos(math.pi * x[1]) + (x[0] * x[1]) ** 3
    elif nr == 7:
        value = 2 * math.pi ** 2 * math.cos(math.pi * x[0]) \
            * math.cos(math.pi * x[1]) \
            - 6 * (x[0] * x[1] ** 3 + x[0] ** 3 * x[1])
    elif nr == 8:
        value = math.pi * math.sin(math.pi * x[0]) \
            * math.cos(math.pi * x[1]) - 3 * x[0] ** 2 * x[1] ** 3
    else:
        raise ValueError("Invalid value for 'nr'")

    return value
