import numpy as np
from scipy.sparse import eye, lil_matrix


def apply_essential(A, f, uess, iess):
    """
    Add effect of essential boundary conditions to right-hand side.

    NOTE: This function modifies the input arguments `A` and `f` in place.

    Parameters
    ----------
    Ain : scipy.sparse.lil_matrix
        System matrix.
    fin : numpy.ndarray
        Right-hand side vector.
    uess : numpy.ndarray
        Vector containing the values for essential boundary conditions.
    iess : numpy.ndarray
        Index of defined essential degrees.

    Returns
    -------
    Aup : scipy.sparse.lil_matrix
        Matrix partition of unknown rows and prescribed columns.
    """

    assert ( uess is not None )

    # initialize some parameters
    nd = f.shape[0]
    npp = iess.shape[0]
    nu = nd - npp

    # find iu
    tmp = np.zeros(nd, dtype=int)
    tmp[iess] = 1
    iu = np.where(tmp == 0)[0]

    # extract Aup
    Aup = A[iu, :][:, iess]

    # modify A
    A[iu[:, None], iess] = lil_matrix((nu, npp))
    A[iess[:, None], iu] = lil_matrix((npp, nu))
    A[iess[:, None], iess] = eye(npp)

    # modify f
    f[iu] -= Aup.dot(uess[iess])
    f[iess] = uess[iess]

    return Aup
