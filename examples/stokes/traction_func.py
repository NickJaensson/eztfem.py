'''
Traction functions for the Stokes example problems.
'''
import numpy as np


def traction_func(nr, x):
    """
    Evaluate a vector function as a function of spatial coordinates.

    Parameters
    ----------
    nr : int
        Function selector.
    x : array_like of shape (ndim,)
        Spatial coordinates, where ndim is the spatial dimension.

    Returns
    -------
    value : numpy.ndarray of shape (2,)
        Vector function value at x.
    """

    value = np.zeros(2)

    if nr == 1:
        value[0] = 1.0
        value[1] = 0.0
    elif nr == 2:
        value[0] = 0.0
        value[1] = x[1]
    else:
        raise ValueError("Invalid value for 'nr'")

    return value
