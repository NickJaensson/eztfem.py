'''
Module for getting the positions of degrees of freedom in nodes.
'''
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
        The requested sequence order of the degrees of freedom in the nodes.
        'ND' : the most inner loop is over the degrees of freedom.
        'DN' : the most inner loop is over the nodal points.
        Default: 'DN'.
        NOTE: the outside loop is always given by the physical quantities.

    Returns
    -------
    pos : list of arrays
        List of arrays containing the positions of the degrees of freedom of
        each physical quantity.
    ndof : list of int
        List of the number of degrees of freedom in `pos` of each physical
        quantity.

    Note
    ----
    In eztfem, the ordering of the full system vector and vectors is always
    NPD. By using the order argument, pos_array returns the indices in the
    requested order. Note, that the default is 'DN', since that is the default
    in build_system.

    Examples
    --------
    >>> pos, ndof = pos_array(problem, nodes, physq=[1, 2], order='ND')

    """

    # Set default optional arguments
    physq = np.arange(problem.nphysq)
    order = 'DN'

    # Override optional arguments if provided
    physq = kwargs.get('physq', physq)
    order = kwargs.get('order', order)

    # Convert physq to a list if an int is supplied
    if isinstance(physq, (int, np.integer)):
        physq = [physq]

    # Convert nodes to a numpy array if an int is supplied
    if isinstance(nodes, (int, np.integer)):
        nodes = np.array([nodes])

    # Convert nodes to a numpy array if a list is supplied
    if isinstance(nodes, list):
        nodes = np.array(nodes)

    # Validate order parameter
    if order not in ['ND', 'DN']:
        raise ValueError(f"Invalid order '{order}'. Must be 'ND' or 'DN'.")

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
        The requested sequence order of the degrees of freedom in the nodes.
        'ND' : the most inner loop is over the degrees of freedom.
        'DN' : the most inner loop is over the nodal points.
        Default: 'DN'.
        NOTE: the outside loop is always given by the physical quantities.

    Returns
    -------
    pos : list of arrays
        List of arrays containing the positions of the degrees of freedom of
        each vector.
    ndof : list of int
        List of the number of degrees of freedom in `pos` of each vector.

    Note
    ----
    In eztfem, the ordering of the full system vector and vectors is always
    NPD. By using the order argument, pos_array returns the indices in the
    requested order. Note, that the default is 'DN', since that is the default
    in build_system.

    Examples
    --------
    >>> pos, ndof = pos_array_vec(problem, nodes, vec=[1, 2], order='ND')

    """

    # Set default optional arguments
    vec = np.arange(problem.nvec)
    order = 'DN'

    # Override optional arguments if provided
    vec = kwargs.get('vec', vec)
    order = kwargs.get('order', order)

    # Convert vec to a numpy array if an int is supplied
    if isinstance(vec, (int, np.integer)):
        vec = np.array([vec])

    # Convert nodes to a numpy array if an int is supplied
    if isinstance(nodes, (int, np.integer)):
        nodes = np.array([nodes])

    # Validate order parameter
    if order not in ['ND', 'DN']:
        raise ValueError(f"Invalid order '{order}'. Must be 'ND' or 'DN'.")

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
