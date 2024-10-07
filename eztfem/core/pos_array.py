import numpy as np


def pos_array(problem, nodes, **kwargs):
    """
    Get the index of the system degrees of freedom in the given nodes.

    Parameters
    ----------
    problem : Problem
        Problem object representing the problem structure.
    nodes : array_like
        An array of node numbers.

    Keyword argument
    ----------------
    physq : array_like, optional
        Array of physical quantity numbers. Default is all physical
        quantities.
    order : str, optional
        The sequence order of the degrees of freedom in the nodes.

    Returns
    -------
    pos : list of arrays
        List of arrays containing the positions of the degrees of freedom of
        each physical quantity.
    ndof : list of int
        List of the number of degrees of freedom in `pos` of each physical
        quantity.

    Examples
    --------
    >>> pos, ndof = get_dof_indices(problem, nodes, physq=[1, 2], order='ND')

    """

    # Set default optional arguments
    physq = np.arange(problem.nphysq)
    order = 'DN'

    # Override optional arguments if provided
    if 'physq' in kwargs:
        physq = kwargs['physq']
    if 'order' in kwargs:
        order = kwargs['order']

    # Convert physq to a list if an int is supplied
    if isinstance(physq, (int, np.integer)):
        physq = [physq]

    # Convert nodes to a numpy array if an int is supplied
    if isinstance(nodes, (int, np.integer)):
        nodes = np.array([nodes])

    pos = [None] * len(physq)
    ndof = np.zeros(len(physq), dtype=int)
    lpos = np.zeros(problem.maxnoddegfd * nodes.shape[0], dtype=int)
    if order == 'ND':
        for i, phq in enumerate(physq):
            dof = 0
            for nodenr in nodes:
                bp = problem.nodnumdegfd[nodenr] + \
                     sum(problem.vec_nodnumdegfd[nodenr+1, :phq]
                         - problem.vec_nodnumdegfd[nodenr, :phq])
                nndof = problem.vec_nodnumdegfd[nodenr+1, phq] \
                    - problem.vec_nodnumdegfd[nodenr, phq]
                lpos[dof:dof+nndof] = np.arange(bp, bp+nndof, dtype=int)
                dof += nndof
            # convert to list: apparently mixing np arrays and lists goes wrong
            pos[i] = lpos[:dof].tolist()
            ndof[i] = dof
    elif order == 'DN':
        for i, phq in enumerate(physq):
            maxdeg = max(problem.vec_nodnumdegfd[nodes+1, phq]
                         - problem.vec_nodnumdegfd[nodes, phq])
            dof = 0
            for deg in range(maxdeg):
                for nodenr in nodes:
                    bp = problem.nodnumdegfd[nodenr] + \
                         sum(problem.vec_nodnumdegfd[nodenr+1, :phq]
                             - problem.vec_nodnumdegfd[nodenr, :phq])
                    nndof = problem.vec_nodnumdegfd[nodenr+1, phq] \
                        - problem.vec_nodnumdegfd[nodenr, phq]
                    if deg < nndof:
                        lpos[dof] = bp + deg
                        dof += 1
            # convert to list: apparently mixing np arrays and lists goes wrong
            pos[i] = lpos[:dof].tolist()
            ndof[i] = dof

    return pos, ndof


def pos_array_vec(problem, nodes, **kwargs):
    """
    Get the index of the degrees of freedom of one or more vectors of special
    structure in the given nodes.

    Parameters
    ----------
    problem : Problem
        Problem object representing the problem structure.
    nodes : array_like
        An array of node numbers.

    Keyword argument
    ----------------
    vec : array_like, optional
        Array of vector numbers. Default is all vectors.
    order : str, optional
        The sequence order of the degrees of freedom in the nodes.

    Returns
    -------
    pos : list of arrays
        List of arrays containing the positions of the degrees of freedom of
        each vector.
    ndof : list of int
        List of the number of degrees of freedom in `pos` of each vector.

    Examples
    --------
    >>> pos, ndof = get_vector_dof_indices(problem, nodes, vec=[1, 2],
                                           order='ND')

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
    ndof = np.zeros(vec.shape[0], dtype=int)
    lpos = np.zeros(problem.maxvecnoddegfd * nodes.shape[0], dtype=int)

    if order == 'ND':
        for i, vc in enumerate(vec):
            dof = 0
            for nodenr in nodes:
                bp = problem.vec_nodnumdegfd[nodenr, vc]
                nndof = problem.vec_nodnumdegfd[nodenr+1, vc] - bp
                lpos[dof:dof+nndof] = np.arange(bp, bp+nndof)
                dof += nndof
            # convert to list: apparently mixing np arrays and lists goes wrong
            pos[i] = lpos[:dof].tolist()
            ndof[i] = dof
    elif order == 'DN':
        for i, vc in enumerate(vec):
            maxdeg = max(problem.vec_nodnumdegfd[nodes+1, vc]
                         - problem.vec_nodnumdegfd[nodes, vc])
            dof = 0
            for deg in range(maxdeg):
                for nodenr in nodes:
                    bp = problem.vec_nodnumdegfd[nodenr, vc]
                    nndof = problem.vec_nodnumdegfd[nodenr+1, vc] \
                        - problem.vec_nodnumdegfd[nodenr, vc]
                    if deg < nndof:
                        lpos[dof] = bp + deg
                        dof += 1
            # convert to list: apparently mixing np arrays and lists goes wrong
            pos[i] = lpos[:dof].tolist()
            ndof[i] = dof

    return pos, ndof
