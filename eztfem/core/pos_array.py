"""Module for getting the positions of degrees of freedom in nodes."""
from typing import Union, List, Tuple
import numpy as np
from numpy.typing import ArrayLike


def pos_array(
    problem,
    nodes: Union[int, ArrayLike],
    physq: Union[int, ArrayLike, None] = None,
    order: str = 'DN'
) -> Tuple[List[List[int]], np.ndarray]:
    """
    Get the index of the system degrees of freedom in the given nodes.

    Parameters
    ----------
    problem : Problem
        Problem object representing the problem structure.
    nodes : array_like
        An array of node numbers.
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
    ndof : ndarray of int
        Array of the number of degrees of freedom in `pos` of each physical
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
    # Validate order parameter
    if order not in ('ND', 'DN'):
        raise ValueError(f"Invalid order '{order}'. Must be 'ND' or 'DN'.")

    # Normalize inputs
    physq_arr = _normalize_physq(physq, problem.nphysq)
    nodes_arr = np.atleast_1d(np.asarray(nodes))

    # Preallocate arrays
    pos = [None] * len(physq_arr)
    ndof = np.zeros(len(physq_arr), dtype=int)
    lpos = np.zeros(problem.maxnoddegfd * len(nodes_arr), dtype=int)

    if order == 'ND':
        _compute_nd_order(problem, nodes_arr, physq_arr, pos, ndof, lpos)
    else:  # order == 'DN'
        _compute_dn_order(problem, nodes_arr, physq_arr, pos, ndof, lpos)

    return pos, ndof


def _normalize_physq(
    physq: Union[int, ArrayLike, None],
    nphysq: int
) -> np.ndarray:
    """Normalize physq parameter to numpy array."""
    if physq is None:
        return np.arange(nphysq)
    return np.atleast_1d(np.asarray(physq))


def _compute_nd_order(
    problem,
    nodes_arr: np.ndarray,
    physq_arr: np.ndarray,
    pos: List,
    ndof: np.ndarray,
    lpos: np.ndarray
) -> None:
    """Compute positions for ND order (degrees of freedom inner loop)."""
    for i, phq in enumerate(physq_arr):
        dof = 0
        for nodenr in nodes_arr:
            base_pos = problem.nodnumdegfd[nodenr] + sum(
                problem.vec_nodnumdegfd[nodenr + 1, :phq] -
                problem.vec_nodnumdegfd[nodenr, :phq]
            )
            num_dof = (
                problem.vec_nodnumdegfd[nodenr + 1, phq] -
                problem.vec_nodnumdegfd[nodenr, phq]
            )
            lpos[dof:dof + num_dof] = np.arange(
                base_pos, base_pos + num_dof, dtype=int
            )
            dof += num_dof
        pos[i] = lpos[:dof].tolist()
        ndof[i] = dof


def _compute_dn_order(
    problem,
    nodes_arr: np.ndarray,
    physq_arr: np.ndarray,
    pos: List,
    ndof: np.ndarray,
    lpos: np.ndarray
) -> None:
    """Compute positions for DN order (nodal points inner loop)."""
    for i, phq in enumerate(physq_arr):
        max_degrees = max(
            problem.vec_nodnumdegfd[nodes_arr + 1, phq] -
            problem.vec_nodnumdegfd[nodes_arr, phq]
        )
        dof = 0
        for deg in range(max_degrees):
            for nodenr in nodes_arr:
                base_pos = problem.nodnumdegfd[nodenr] + sum(
                    problem.vec_nodnumdegfd[nodenr + 1, :phq] -
                    problem.vec_nodnumdegfd[nodenr, :phq]
                )
                num_dof = (
                    problem.vec_nodnumdegfd[nodenr + 1, phq] -
                    problem.vec_nodnumdegfd[nodenr, phq]
                )
                if deg < num_dof:
                    lpos[dof] = base_pos + deg
                    dof += 1
        pos[i] = lpos[:dof].tolist()
        ndof[i] = dof


def pos_array_vec(
    problem,
    nodes: Union[int, ArrayLike],
    vec: Union[int, ArrayLike, None] = None,
    order: str = 'DN'
) -> Tuple[List[List[int]], np.ndarray]:
    """
    Get the index of the degrees of freedom of one or more vectors of special
    structure in the given nodes.

    Parameters
    ----------
    problem : Problem
        Problem object representing the problem structure.
    nodes : array_like
        An array of node numbers.
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
    ndof : ndarray of int
        Array of the number of degrees of freedom in `pos` of each vector.

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
    # Validate order parameter
    if order not in ('ND', 'DN'):
        raise ValueError(f"Invalid order '{order}'. Must be 'ND' or 'DN'.")

    # Normalize inputs
    vec_arr = _normalize_vec(vec, problem.nvec)
    nodes_arr = np.atleast_1d(np.asarray(nodes))

    # Preallocate arrays
    pos = [None] * len(vec_arr)
    ndof = np.zeros(len(vec_arr), dtype=int)
    lpos = np.zeros(problem.maxvecnoddegfd * len(nodes_arr), dtype=int)

    if order == 'ND':
        _compute_vec_nd_order(problem, nodes_arr, vec_arr, pos, ndof, lpos)
    else:  # order == 'DN'
        _compute_vec_dn_order(problem, nodes_arr, vec_arr, pos, ndof, lpos)

    return pos, ndof


def _normalize_vec(
    vec: Union[int, ArrayLike, None],
    nvec: int
) -> np.ndarray:
    """Normalize vec parameter to numpy array."""
    if vec is None:
        return np.arange(nvec)
    return np.atleast_1d(np.asarray(vec))


def _compute_vec_nd_order(
    problem,
    nodes_arr: np.ndarray,
    vec_arr: np.ndarray,
    pos: List,
    ndof: np.ndarray,
    lpos: np.ndarray
) -> None:
    """Compute positions for ND order (degrees of freedom inner loop)."""
    for i, vc in enumerate(vec_arr):
        dof = 0
        for nodenr in nodes_arr:
            base_pos = problem.vec_nodnumdegfd[nodenr, vc]
            num_dof = problem.vec_nodnumdegfd[nodenr + 1, vc] - base_pos
            lpos[dof:dof + num_dof] = np.arange(base_pos, base_pos + num_dof)
            dof += num_dof
        pos[i] = lpos[:dof].tolist()
        ndof[i] = dof


def _compute_vec_dn_order(
    problem,
    nodes_arr: np.ndarray,
    vec_arr: np.ndarray,
    pos: List,
    ndof: np.ndarray,
    lpos: np.ndarray
) -> None:
    """Compute positions for DN order (nodal points inner loop)."""
    for i, vc in enumerate(vec_arr):
        max_degrees = max(
            problem.vec_nodnumdegfd[nodes_arr + 1, vc] -
            problem.vec_nodnumdegfd[nodes_arr, vc]
        )
        dof = 0
        for deg in range(max_degrees):
            for nodenr in nodes_arr:
                base_pos = problem.vec_nodnumdegfd[nodenr, vc]
                num_dof = (
                    problem.vec_nodnumdegfd[nodenr + 1, vc] -
                    problem.vec_nodnumdegfd[nodenr, vc]
                )
                if deg < num_dof:
                    lpos[dof] = base_pos + deg
                    dof += 1
        pos[i] = lpos[:dof].tolist()
        ndof[i] = dof
