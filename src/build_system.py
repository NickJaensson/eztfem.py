import numpy as np
from scipy.sparse import lil_matrix
from src.pos_array import pos_array
from src.pos_array_vec import pos_array_vec

def build_system(mesh, problem, element, user, **kwargs):
    """
    Build the system matrix

    Parameters:
    mesh: Mesh structure
    problem: Problem object
    element: Function handle to the element function routine
    user: Can be used by the user for transferring data to the element routine
    **kwargs: Optional arguments
              - physqrow: array of physical quantity numbers for the rows of the matrix and for the right-hand side vector
              - physqcol: array of physical quantity numbers for the columns of the matrix
              - order: the sequence order of the degrees of freedom on element level
              - posvectors: supply the position of vectors to the element routine
    
    Returns:
    A: System matrix
    f: System vector
    """

    # Set default optional arguments
    physqrow = np.arange(problem.nphysq,dtype=int)
    physqcol = np.arange(problem.nphysq,dtype=int)
    order = 'DN'
    posvectors = False

    # Override optional arguments if provided
    if 'physqrow' in kwargs:
        physqrow = kwargs['physqrow']
    if 'physqcol' in kwargs:
        physqcol = kwargs['physqcol']
    if 'order' in kwargs:
        order = kwargs['order']
    if 'posvectors' in kwargs:
        posvectors = kwargs['posvectors']

    rowcolequal = np.array_equal(physqrow, physqcol)

    n = problem.numdegfd
    A = lil_matrix((n, n))
    f = np.zeros(n)

    # Start assembly loop over elements
    for elem in range(mesh.nelem):

        posrow, _ = pos_array(problem, mesh.topology[:,elem].T, order=order)
 
        posr = np.hstack([posrow[i] for i in physqrow]) # indexing a list using another list

        if rowcolequal:
            posc = posr
        else:
            poscol, _ = pos_array(problem, mesh.topology[:,elem].T, physq=physqcol, order=order)
            posc = np.hstack(poscol)

        coor = mesh.coor[mesh.topology[:,elem],:]

        if posvectors:
            posvec, _ = pos_array_vec(problem, mesh.topology[:,elem].T, order=order)
            elemmat, elemvec = element(elem, coor, user, posrow, posvec)
        else:
            elemmat, elemvec = element(elem, coor, user, posrow)

        # [:,None] needed for proper broadcasting by adding an axis of dim 1
        A[posr[:,None], posc] += elemmat
        f[posr] += elemvec

    return A, f
