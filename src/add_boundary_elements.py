import numpy as np
from src.pos_array import pos_array
from src.pos_array_vec import pos_array_vec

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
    curve = kwargs.get('curve', 0)
    Ap = kwargs.get('Ap', None)
    physqrow = np.arange(1, problem.nphysq + 1,dtype=int)
    physqcol = np.arange(1, problem.nphysq + 1,dtype=int)
    order = kwargs.get('order', 'DN')
    posvectors = kwargs.get('posvectors', 0)

    if curve == 0:
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
        # positions in global system (Replace 'pos_array' with appropriate Python function)
        posrow = pos_array(problem, mesh.curves[curve].topology[elem, 2, :], order=order)
        print(posrow)
        posr = np.concatenate(posrow[physqrow])

        if mat:
            if rowcolequal:
                posc = posr
            else:
                poscol = pos_array(problem, mesh.curves[curve].topology[elem, 2, :], physq=physqcol, order=order)
                posc = np.concatenate(poscol)

        coor = mesh.coor[mesh.curves[curve].topology[elem, 2, :], :]

        if posvectors:
            posvec = pos_array_vec(problem, mesh.curves[curve].topology[elem, 2, :], order=order)
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
