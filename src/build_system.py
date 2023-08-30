import numpy as np
from scipy.sparse import lil_matrix
from src.pos_array import pos_array
from src.pos_array_vec import pos_array_vec

def build_system(mesh, problem, element, user, **kwargs):
    """
    BUILD_SYSTEM  Build the system matrix
      [ A, f ] = BUILD_SYSTEM ( mesh, problem, element, user, 'option1', value1, .... )
      input:
        mesh: mesh structure
        mesh: problem structure
        element: function handle to the element function routine
        user: can be used by the user for transferring data to the element routine
      optional arguments:
        string, value couples to set options:
        'physqrow' array of physical quantity numbers for the rows of the matrix
                   and for the right-hand side vector.
                   default: all physical quantities
        'physqcol' array of physical quantity numbers for the columns of the matrix
                   default: all physical quantities
        'order'  the sequence order of the degrees of freedom on element level:
                 'ND' : the most inner loop is over the degrees of freedom
                 'DN' : the most inner loop is over the nodal points
                 default = 'DN'  
                 NOTE: the outside loop is always given by the physical quantities.
        'posvectors' supply the position of vectors to the element routine
                 default=0
        For example:
          [ A, f ] = build_system ( mesh, problem, @element, user, 'order', 'ND' )
        to change the order.
      output:
        A: the system matrix. 
        f: the system vector. 
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
