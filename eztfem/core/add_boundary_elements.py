import numpy as np
from .pos_array import pos_array
from .pos_array_vec import pos_array_vec


def add_boundary_elements(mesh, problem, f, element, user, curve, **kwargs):
    """
    Add boundary elements to the system vector and optionally to the system
    matrix.

    NOTE: This function modifies the input arguments `f` and optionally `A` in
    place.

    Parameters
    ----------
    mesh : Mesh
        Mesh object.
    problem : Problem
        Problem object.
    f : numpy.ndarray
        "Previous" system vector input. Result added to this vector.
    element : callable
        Function handle to the element function routine.
    user : Any
        Can be used by the user for transferring data to the element routine.
    curve : int
        Build on given curve number.

    Keyword arguments
    -----------------
    A : numpy.ndarray
        "Previous" system matrix. If present, element matrices will be
        needed and added to this matrix.
    physqrow : numpy.ndarray
        Array of physical quantity numbers for the rows of the matrix
        and for the right-hand side vector. Default: All physical
        quantities.
    physqcol : numpy.ndarray
        Array of physical quantity numbers for the columns of the matrix.
        Default: All physical quantities.
    order : str
        The sequence order of the degrees of freedom on element level:
        'ND' or 'DN'. Default = 'DN'
    posvectors : bool
        Supply the position of vectors to the element routine.
        Default=False

    """

    # Default optional arguments
    A = kwargs.get('A', None)
    physqrow = kwargs.get('physqrow', np.arange(problem.nphysq, dtype=int))
    physqcol = kwargs.get('physqcol', np.arange(problem.nphysq, dtype=int))
    order = kwargs.get('order', 'DN')
    posvectors = kwargs.get('posvectors', False)

    mat_present = A is not None

    rowcolequal = np.array_equal(physqrow, physqcol)

    # Assemble loop over elements
    for elem in range(mesh.curves[curve].nelem):

        nodes = mesh.curves[curve].topology[:, elem, 1].T

        posrow, _ = pos_array(problem, nodes, order=order)

        # indexing a list using another list
        posr = np.hstack([posrow[i] for i in physqrow])

        if mat_present:
            if rowcolequal:
                posc = posr
            else:
                poscol, _ = pos_array(problem, nodes, physq=physqcol,
                                      order=order)
                posc = np.concatenate(poscol)

        coor = mesh.coor[mesh.curves[curve].topology[:, elem, 1], :]

        # NEED FURTHER CHECKIING IN FULL CODE IF POS_ARRAY AND POS_ARRAY_VEC
        # OUTPUT IS ALWAYS DONE WITH A _

        if posvectors:
            posvec, _ = pos_array_vec(problem, nodes, order=order)
            if mat_present:
                elemvec, elemmat = element(elem, coor, user, posrow, posvec)
                A[posr[:, None], posc] += elemmat
            else:
                elemvec = element(elem, coor, user, posrow, posvec)
        else:
            if mat_present:
                elemvec, elemmat = element(elem, coor, user, posrow)
                A[posr[:, None], posc] += elemmat
            else:
                elemvec = element(elem, coor, user, posrow)

        f[posr] += elemvec
