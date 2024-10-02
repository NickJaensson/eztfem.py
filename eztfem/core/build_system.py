import numpy as np
from scipy.sparse import lil_matrix
from .pos_array import pos_array
from .pos_array_vec import pos_array_vec


def build_system(mesh, problem, element, user, **kwargs):
    """
    Build the system matrix and right hand side.

    Parameters
    ----------
    mesh : Mesh
        Mesh object.
    problem : Problem
        Problem object.
    element : callable
        Function handle to the element function routine.
    user : User
        User object to pass parameters and data to the element routine.

    Keyword arguments
    -----------------
    physqrow : numpy.ndarray, optional
        Array of physical quantity numbers for the rows of the matrix
        and for the right-hand side vector. Default: all physical
        quantities.
    physqcol : numpy.ndarray, optional
        Array of physical quantity numbers for the columns of the matrix.
        Default: all physical quantities.
    order : str, optional
        The sequence order of the degrees of freedom on element level:
        'ND' : the most inner loop is over the degrees of freedom.
        'DN' : the most inner loop is over the nodal points.
        Default: 'DN'.
        NOTE: the outside loop is always given by the physical quantities.
    posvectors : bool, optional
        Supply the position of vectors to the element routine.
        Default: False.

    Returns
    -------
    A : numpy.ndarray
        The system matrix.
    f : numpy.ndarray
        The right hand side.

    Examples
    --------
    To change the order:
    >>> A, f = build_system(mesh, problem, element, user, order='ND')

    """

    # Set default optional arguments
    physqrow = kwargs.get('physqrow', np.arange(problem.nphysq, dtype=int))
    physqcol = kwargs.get('physqcol', np.arange(problem.nphysq, dtype=int))
    order = kwargs.get('order', 'DN')
    posvectors = kwargs.get('posvectors', False)

    rowcolequal = np.array_equal(physqrow, physqcol)

    n = problem.numdegfd
    A = lil_matrix((n, n))
    f = np.zeros(n)

    # Start assembly loop over elements
    for elem in range(mesh.nelem):

        nodes = mesh.topology[:, elem].T

        posrow, _ = pos_array(problem, nodes, order=order)

        # indexing a list using another list
        posr = np.hstack([posrow[i] for i in physqrow])

        if rowcolequal:
            posc = posr
        else:
            poscol, _ = pos_array(problem, nodes, physq=physqcol, order=order)
            posc = np.hstack(poscol)

        coor = mesh.coor[mesh.topology[:, elem], :]

        if posvectors:
            posvec, _ = pos_array_vec(problem, nodes, order=order)
            elemmat, elemvec = element(elem, coor, user, posrow, posvec)
        else:
            elemmat, elemvec = element(elem, coor, user, posrow)

        # [:,None] needed for proper broadcasting by adding an axis of dim 1
        A[posr[:, None], posc] += elemmat
        f[posr] += elemvec

    return A, f
