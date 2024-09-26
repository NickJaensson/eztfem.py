import numpy as np

def pos_array_vec(problem, nodes, **kwargs):
    """
    Get the index of the degrees of freedom of one or more vectors of special structure in the given nodes.

    Parameters
    ----------
    problem : Problem
        Problem object representing the problem structure.
    nodes : array_like
        An array of node numbers.
    **kwargs : optional
        Additional options.
        - vec : array_like, optional
            Array of vector numbers. Default is all vectors.
        - order : str, optional
            The sequence order of the degrees of freedom in the nodes.

    Returns
    -------
    pos : list of arrays
        List of arrays containing the positions of the degrees of freedom of each vector.
    ndof : list of int
        List of the number of degrees of freedom in `pos` of each vector.

    Examples
    --------
    >>> pos, ndof = get_vector_dof_indices(problem, nodes, vec=[1, 2], order='ND')
    """

    # Set default optional arguments
    vec = np.arange(problem.nvec)
    order = 'DN'
    
    # Override optional arguments if provided
    if 'vec' in kwargs:
        vec = kwargs['vec']
    if 'order' in kwargs:
        order = kwargs['order']

    # Convert vec to a numpy array if an int is supplied
    if isinstance(vec, (int, np.integer)):
        vec = np.array([vec])

    # Convert nodes to a numpy array if an int is supplied
    if isinstance(nodes, (int, np.integer)):
        nodes = np.array([nodes])

    pos = [None] * vec.shape[0]
    ndof = np.zeros(vec.shape[0],dtype=int)
    lpos = np.zeros(problem.maxvecnoddegfd * nodes.shape[0],dtype=int)

    if order == 'ND':
        for i, vc in enumerate(vec):
            dof = 0
            for nodenr in nodes:
                bp = problem.vec_nodnumdegfd[nodenr, vc]
                nndof = problem.vec_nodnumdegfd[nodenr+1, vc] - bp
                lpos[dof:dof+nndof] = np.arange(bp, bp+nndof)
                dof += nndof
            pos[i] = lpos[:dof].tolist() # convert to list: apparently mixing np arrays and lists goes wrong!
            ndof[i] = dof
    elif order == 'DN':
        for i, vc in enumerate(vec):
            maxdeg = max(problem.vec_nodnumdegfd[nodes+1, vc] - problem.vec_nodnumdegfd[nodes, vc])
            dof = 0
            for deg in range(maxdeg):
                for nodenr in nodes:
                    bp = problem.vec_nodnumdegfd[nodenr, vc]
                    nndof = problem.vec_nodnumdegfd[nodenr+1, vc] - problem.vec_nodnumdegfd[nodenr, vc]
                    if deg < nndof:
                        lpos[dof] = bp + deg
                        dof += 1
            pos[i] = lpos[:dof].tolist() # convert to list: apparently mixing np arrays and lists goes wrong!
            ndof[i] = dof

    return pos, ndof 
