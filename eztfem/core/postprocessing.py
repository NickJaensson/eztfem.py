"""Module for postprocessing routines"""
from typing import Any, Callable

from .meshgen import Mesh
from .problem import Problem
from .user import User
from .pos_array import pos_array, pos_array_vec


def integrate_boundary_elements(
    mesh: Mesh,
    problem: Problem,
    element: Callable,
    user: User,
    **kwargs: Any,
) -> float:
    """
    Integrate boundary elements.

    Parameters
    ----------
    mesh : Mesh
        Mesh object.
    problem : Problem
        Problem object.
    element : Callable
        Function handle to the element function routine.
    user : User
        User object.

    Keyword arguments
    -----------------
    curve : int
        Build on given curve number.
    order : {'ND', 'DN'}, default='DN'
        The sequence order of the degrees of freedom on element level:
        - 'ND' : the most inner loop is over the degrees of freedom.
        - 'DN' : the most inner loop is over the nodal points.
        NOTE: The outside loop is always given by the physical quantities.
    posvectors : bool, default=False
        Supply the position of vectors to the element routine.

    Returns
    -------
    resultsum : float
        The sum over all elements of the integral on elements.

    Note
    ----
    In eztfem, the ordering of the full system vector and vectors is always
    NPD. With the order argument, we can specify the ordering on how the
    element routine is written.

    Examples
    --------
    >>> resultsum = integrate_boundary_elements(mesh, problem, element, user)

    """
    curve = kwargs.get("curve", 0)
    order = kwargs.get("order", "DN")
    posvectors = kwargs.get("posvectors", False)

    if curve == 0:
        raise ValueError("Argument 'curve' is required and must be non-zero.")

    boundary_curve = mesh.curves[curve]
    resultsum = 0

    for elem in range(boundary_curve.nelem):
        elem_topology = boundary_curve.topology[:, elem, 1]
        nodes = elem_topology.T

        # Positions in global system
        pos, _ = pos_array(problem, nodes, order=order)
        coor = mesh.coor[elem_topology, :]

        # Call element routine with appropriate arguments
        element_args = (elem, coor, user, pos)
        if posvectors:
            posvec, _ = pos_array_vec(problem, nodes, order=order)
            resultsum += element(*element_args, posvec)
        else:
            resultsum += element(*element_args)

    return resultsum
