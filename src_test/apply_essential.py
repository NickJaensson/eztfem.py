import numpy as np
from scipy.sparse import eye, lil_matrix

def apply_essential(Ain, fin, uess, iess):
    """
    Apply essential degrees of freedom.

    Parameters:
        Ain: system matrix
        fin: right-hand side
        uess: a system vector containing the values for essential bc
        iess: index of defined essential degrees

    Returns:
        A, f: the system matrix and right-hand side with essential bc applied
        Aup: matrix partition of unknown rows and prescribed columns.
    """
    
    # Copy the system matrix and vector
    A = Ain.copy()
    f = np.array(fin)
    
    nd = fin.shape[0]
    npp = iess.shape[0]
    nu = nd - npp

    # Find iu
    tmp = np.zeros(nd, dtype=int)
    tmp[iess] = 1
    iu = np.where(tmp == 0)[0]

    # Extract Aup

    Aup = Ain[iu,:][:,iess]

    # Modify A
    A[iu[:,None], iess] = lil_matrix((nu, npp))
    A[iess[:,None], iu] = lil_matrix((npp, nu))
    A[iess[:,None], iess] = eye(npp)

    # Modify f
    f[iu] -= Aup.dot(uess[iess])
    f[iess] = uess[iess]

    return A, f, Aup
