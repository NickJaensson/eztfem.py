import typing
import numpy as np
from .tp import ArrayLike, Order

if typing.TYPE_CHECKING:
    from .problem import Problem


def pos_array(problem: "Problem", nodes: ArrayLike[int, np.integer], *,
              physq: ArrayLike[int, np.integer] | None = None,
              order: Order = "DN"):
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

    if physq is None:
        physq = np.arange(problem.nphysq)

    # Convert physq to a list if an int is supplied
    if isinstance(physq, (int, np.integer)):
        physq = typing.cast(list[int], [physq])

    # Convert nodes to a numpy array if an int is supplied
    if isinstance(nodes, (int, np.integer)):
        nodes = np.array([nodes])

    # Convert nodes to a numpy array if a list is supplied
    elif isinstance(nodes, typing.Sequence):
        nodes = np.array(nodes)

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

    # NOTE: Typechecker infers pos as list[None] because of the preallocation.
    return typing.cast(list[list[int]], pos), ndof


def pos_array_vec(problem: "Problem", nodes: ArrayLike[int, np.integer], *,
                  vec: ArrayLike[int, np.integer] | None = None,
                  order: Order = "DN"):
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

    if vec is None:
        vec = np.arange(problem.nvec)

    # Convert vec to a numpy array if an int is supplied
    if isinstance(vec, (int, np.integer)):
        vec = np.array([vec])

    # Convert vec to a numpy array if a list is supplied
    elif isinstance(vec, typing.Sequence):
        vec = np.array(vec)

    # Convert nodes to a numpy array if an int is supplied
    if isinstance(nodes, (int, np.integer)):
        nodes = np.array([nodes])

    # Convert nodes to a numpy array if a list is supplied
    elif isinstance(nodes, typing.Sequence):
        nodes = np.array(nodes)

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

    # NOTE: Typechecker infers pos as list[None] because of the preallocation.
    return typing.cast(list[list[int]], pos), ndof
