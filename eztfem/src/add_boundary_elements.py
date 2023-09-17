import numpy as np
from eztfem.src.pos_array import pos_array
from eztfem.src.pos_array_vec import pos_array_vec

def add_boundary_elements(mesh, problem, fp, element, user, **kwargs):
    """
    Add boundary elements to the system vector and optionally to the system matrix.
    
    Parameters:
        mesh : dict
            Mesh structure.
        problem : dict
            Problem structure.
        fp : np.ndarray
            "Previous" system vector input. Result added to this vector.
        element : callable
            Function handle to the element function routine.
        user : Any
            Can be used by the user for transferring data to the element routine.
        kwargs : dict
            Optional arguments:
            curve : int
                Build on given curve number. (Required)
            Ap : np.ndarray
                "Previous" system matrix. If present, element matrices will be needed
                and added to this matrix.
            physqrow : np.ndarray
                Array of physical quantity numbers for the rows of the matrix
                and for the right-hand side vector.
                Default: All physical quantities.
            physqcol : np.ndarray
                Array of physical quantity numbers for the columns of the matrix.
                Default: All physical quantities.
            order : str
                The sequence order of the degrees of freedom on element level:
                'ND' or 'DN'. Default = 'DN'
            posvectors : int
                Supply the position of vectors to the element routine. Default=0

    Returns:
        f : np.ndarray
            The system vector.
        A : np.ndarray
            The system matrix (optionally).
    """

    # Default optional arguments
    curve = kwargs.get('curve', -1)
    Ap = kwargs.get('Ap', None)
    physqrow = kwargs.get('physqrow',np.arange(problem.nphysq,dtype=int))
    physqcol = kwargs.get('physqcol',np.arange(problem.nphysq,dtype=int))
    order = kwargs.get('order', 'DN')
    posvectors = kwargs.get('posvectors', False)

    if curve == -1:
        raise ValueError('Argument curve is missing.')

    # Make a copy of fp to f
    f = fp.copy()
    if Ap is not None:
        A = Ap.copy()
        mat = True
    else:
        mat = False

    rowcolequal = np.array_equal(physqrow, physqcol)

    # Assemble loop over elements
    for elem in range(mesh.curves[curve].nelem):

        posrow, _ = pos_array(problem, mesh.curves[curve].topology[:, elem, 1].T, order=order)

        posr = np.hstack([posrow[i] for i in physqrow]) # indexing a list using another list

        if mat:
            if rowcolequal:
                posc = posr
            else:
                poscol, _ = pos_array(problem, mesh.curves[curve].topology[:, elem, 1].T, physq=physqcol, order=order)
                posc = np.concatenate(poscol)

        coor = mesh.coor[mesh.curves[curve].topology[:, elem, 1], :]

########### NEED FURTHER CHECKIING IN FULL CODE IF POS_ARRAY AND POS_ARRAY_VEC OUTPUT IS ALWAYS DONE WITH A _ !!! #############

        if posvectors:
            posvec, _ = pos_array_vec(problem, mesh.curves[curve].topology[:, elem, 1].T, order=order)
            if mat:
                elemvec, elemmat = element(elem, coor, user, posrow, posvec)
                A[posr[:, None], posc] += elemmat
            else:
                elemvec = element(elem, coor, user, posrow, posvec)
        else:
            if mat:
                elemvec, elemmat = element(elem, coor, user, posrow)
                A[posr[:, None], posc] += elemmat
            else:
                elemvec = element(elem, coor, user, posrow)

        f[posr] += elemvec

    if mat:
        return f, A
    else:
        return f
