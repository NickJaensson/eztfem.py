import typing
from .pos_array import pos_array, pos_array_vec
from .tp import Order, UserT, ElementRoutine, PosvecElementRoutine

if typing.TYPE_CHECKING:
    from .meshgen import Mesh
    from .problem import Problem


@typing.overload
def integrate_boundary_elements(mesh: "Mesh", problem: "Problem",
                                element: PosvecElementRoutine[UserT, float],
                                user: UserT, *, curve: int = 0,
                                order: Order = "DN",
                                posvectors: typing.Literal[True]) -> float:
    ...

@typing.overload
def integrate_boundary_elements(mesh: "Mesh", problem: "Problem",
                                element: ElementRoutine[UserT, float], user: UserT,
                                *, curve: int = 0, order: Order = "DN",
                                posvectors: typing.Literal[False] = False) -> float:
    ...

def integrate_boundary_elements(mesh: "Mesh", problem: "Problem",
                                element: typing.Callable[..., float],
                                user: UserT, *, curve: int = 0,
                                order: Order = "DN",
                                posvectors: bool = False) -> float:
    """
    Integrate boundary elements.

    Parameters
    ----------
    mesh : Mesh
        Mesh object.
    problem : Problem
        Problem object.
    element : function
        Function handle to the element function routine.
    user : any
        Can be used by the user for transferring data to the element routine.

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

    Notes
    -----
    Currently, 'curve' is the only possibility for now.

    Examples
    --------
    >>> resultsum = integrate_boundary_elements(mesh, problem, element, user)

    """

    if curve == 0:
        raise ValueError('Argument curve is missing.')

    # Start summing loop over elements
    resultsum = 0

    for elem in range(mesh.curves[curve].nelem):

        nodes = mesh.curves[curve].topology[:, elem, 1].T

        # Positions in global system
        pos, _ = pos_array(problem, nodes, order=order)
        coor = mesh.coor[mesh.curves[curve].topology[:, elem, 1], :]

        if posvectors:
            posvec, _ = pos_array_vec(problem, nodes, order=order)
            resultsum += element(elem, coor, user, pos, posvec)
        else:
            resultsum += element(elem, coor, user, pos)

    return resultsum
