'''
Supporting functions for the Stokes example problems.
'''
import math


def func(nr, x):
    """
    Evaluate a scalar function as a function of spatial coordinates.

    Parameters
    ----------
    nr : int1
        Function selector.
    x : numpy.ndarray of shape (ndim,)
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
    else:
        raise ValueError("Invalid value for 'nr'")

    return value
