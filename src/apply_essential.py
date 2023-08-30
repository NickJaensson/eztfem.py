import numpy as np
from scipy.sparse import eye, lil_matrix

def apply_essential(Ain, fin, uess, iess):
    """
    APPLY_ESSENTIAL  apply essential degrees of freedom
      [ A, f, Aup ] = APPLY_ESSENTIAL ( Ain, fin, iess )
      input:
        Ain: system matrix
        fin: right-hand side.
        uess: a system vector containing the values for essential bc
        iess: index of defined essential degrees
      output:
        A, f: the system matrix and right-hand side with ess bc applied
        Aup: matrix partition of unknown rows and prescribed columns.
    """
    
    # copy the system matrix and vector
    A = Ain.copy()
    f = np.array(fin)
    
    nd = fin.shape[0]
    npp = iess.shape[0]
    nu = nd - npp

    # find iu
    tmp = np.zeros(nd, dtype=int)
    tmp[iess] = 1
    iu = np.where(tmp == 0)[0]

    # extract Aup
    Aup = Ain[iu,:][:,iess]

    # modify A
    A[iu[:,None], iess] = lil_matrix((nu, npp))
    A[iess[:,None], iu] = lil_matrix((npp, nu))
    A[iess[:,None], iess] = eye(npp)

    # modify f
    f[iu] -= Aup.dot(uess[iess])
    f[iess] = uess[iess]

    return A, f, Aup
