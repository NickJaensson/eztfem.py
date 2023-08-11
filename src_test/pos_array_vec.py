import numpy as np

def pos_array_vec(problem, nodes, **kwargs):
    """
    Get the index of the degrees of freedom of one or more vectors of 
    special structure in the given nodes.

    Parameters:
    problem: object representing the problem structure
    nodes: An array of node numbers
    **kwargs: Optional arguments
        - vec: Array of vector numbers, default: all vectors
        - order: The sequence order of the degrees of freedom in the nodes
    
    Returns:
    pos: List of arrays containing the positions of the degrees of freedom of each vector
    ndof: List of the number of degrees in pos of each vector
    """

    # Set default optional arguments
    vec = np.arange(problem.nvec)
    order = 'DN'
    
    # Override optional arguments if provided
    if 'vec' in kwargs:
        vec = kwargs['vec']
    if 'order' in kwargs:
        order = kwargs['order']

    pos = [None] * len(vec)
    ndof = np.zeros(len(vec))
    lpos = np.zeros(problem.maxvecnoddegfd * len(nodes))

    if order == 'ND':
        for i, vc in enumerate(vec):
            dof = 0
            for nodenr in nodes:
                bp = problem.vec_nodnumdegfd[nodenr, vc]
                nndof = problem.vec_nodnumdegfd[nodenr+1, vc] - bp
                lpos[dof:dof+nndof] = np.arange(bp+1, bp+nndof+1)
                dof += nndof
            pos[i] = lpos[:dof]
            ndof[i] = dof
    elif order == 'DN':
        for i, vc in enumerate(vec):
            maxdeg = max(problem.vec_nodnumdegfd[nodes+1, vc] - problem.vec_nodnumdegfd[nodes, vc])
            dof = 0
            for deg in range(1, maxdeg+1):
                for nodenr in nodes:
                    bp = problem.vec_nodnumdegfd[nodenr, vc]
                    nndof = problem.vec_nodnumdegfd[nodenr+1, vc] - problem.vec_nodnumdegfd[nodenr, vc]
                    if deg <= nndof:
                        lpos[dof] = bp + deg
                        dof += 1
            pos[i] = lpos[:dof] - 1 # Python 0-based indexing
            ndof[i] = dof

    return pos, ndof 
