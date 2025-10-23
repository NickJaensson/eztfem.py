import typing
import numpy as np

from .pos_array import pos_array, pos_array_vec
from .tp import ArrayLike, Order

if typing.TYPE_CHECKING:
    from .meshgen import Mesh
    from .problem import Problem
    from .user import User


class Vector:
    """
    A class used to represent a Vector.

    Attributes
    ----------
    vec : int
        The vector value.
    u : np.ndarray
        An array to store vector data.

    """

    def __init__(self, vec: int = 0):
        """
        Initializes the Vector object with the given attributes. The stored
        vecor data (u) is initialized as an empty numpy array (np.array([])).

        Parameters for initialization
        -----------------------------
        vec : int or float, optional
            The vector value (default is 0).

        """
        self.vec = vec
        self.u = np.array([])

    def __eq__(self, other: "Vector"):  # type: ignore
        """
        Checks equivalence of two Vector objects (overloads == sign).

        Parameters
        ----------
        other : Vector
            The other Vector object to compare with.

        Returns
        -------
        bool
            True if the Vector objects are equivalent, False otherwise.

        Notes
        -----
        See NOTE_ON_COMPARING_ARRAYS.md for the use of np.squeeze.

        """
        check = [self.vec == other.vec,
                 np.allclose(np.squeeze(self.u), other.u, atol=1e-12, rtol=0)]
        return all(check)


@typing.overload
def deriv_vector(mesh: "Mesh", problem: "Problem",
                 element: typing.Callable[[int, np.ndarray, "User", typing.Sequence[typing.Sequence[int]], typing.Sequence[typing.Sequence[int]]], np.ndarray],
                 user: "User", *, vec: int | None = None,
                 order: Order = "DN",
                 posvectors: typing.Literal[True]) -> Vector:
    ...

@typing.overload
def deriv_vector(mesh: "Mesh", problem: "Problem",
                 element: typing.Callable[[int, np.ndarray, "User", typing.Sequence[typing.Sequence[int]]], np.ndarray],
                 user: "User", *, vec: int | None = None,
                 order: Order = "DN",
                 posvectors: typing.Literal[False] = False) -> Vector:
    ...

def deriv_vector(mesh: "Mesh", problem: "Problem",
                 element: typing.Callable[..., typing.Any], user: "User", *,
                 vec: int | None = None,
                 order: Order = "DN",
                 posvectors: bool = False):
    """
    Derive a vector of special structure.

    Parameters
    ----------
    mesh : Mesh
        Mesh structure.
    problem : Problem
        Problem structure.
    element : callable
        Function handle to the element function routine.
    user : Any
        User object to pass parameters and data to the element routine.

    Keyword arguments
    -----------------
    vec : int, optional
        The vector number. Default is problem['nphysq'].
    order : str, optional
        The sequence order of the degrees of freedom on element level.
        'ND' : the most inner loop is over the degrees of freedom.
        'DN' : the most inner loop is over the nodal points.
        Default is 'DN'.
    posvectors : bool, optional
        Supply the position of vectors to the element routine.
        Default: False.

    Returns
    -------
    v : Vector
        Vector object containing the derived data.

    """

    # Set default optional arguments
    # NOTE: default vec is first vec after nphysq (starting at zero)
    if vec is None:
        vec = problem.nphysq

    v = Vector(vec=vec)
    n = problem.vec_numdegfd[vec]
    v.u = np.zeros(n)
    w = np.zeros(n)  # weights

    # Start assemble loop over elements
    for elem in range(mesh.nelem):
        # Positions in the global system
        pos, _ = pos_array(problem, mesh.topology[:, elem].T, order=order)
        posvec, _ = pos_array_vec(problem, mesh.topology[:, elem].T, vec=vec,
                                  order=order)
        posv = np.array(posvec).flatten()
        coor = mesh.coor[mesh.topology[:, elem], :]

        if posvectors:
            posvec, _ = pos_array_vec(problem, mesh.topology[:, elem].T,
                                      order=order)
            elemvec = element(elem, coor, user, pos, posvec)
        else:
            elemvec = element(elem, coor, user, pos)

        v.u[posv] += elemvec
        w[posv] += 1.0

    # Apply averaging by the total weight
    in_ = w > 0
    v.u[in_] = v.u[in_] / w[in_]

    return v


@typing.overload
def fill_system_vector(mesh: "Mesh", problem: "Problem",
                       geometry: typing.Literal["nodes", "points", "curves"],
                       numbers: ArrayLike[int, np.integer],
                       func: typing.Callable[[int, np.ndarray], float], *,
                       funcnr: int = 0, physq: int = 0, degfd: int = 0,
                       f: np.ndarray) -> None:
    ...


@typing.overload
def fill_system_vector(mesh: "Mesh", problem: "Problem",
                       geometry: typing.Literal["nodes", "points", "curves"],
                       numbers: ArrayLike[int, np.integer],
                       func: typing.Callable[[int, np.ndarray], float], *,
                       funcnr: int = 0, physq: int = 0, degfd: int = 0,
                       f: None = None) -> np.ndarray:
    ...


# FIXME: defaults in code are 0; docstring says 1. Also, docstring says func
#        only takes a coordinate vector but it also takes funcnr.
def fill_system_vector(mesh: "Mesh", problem: "Problem",
                       geometry: typing.Literal["nodes", "points", "curves"],
                       numbers: ArrayLike[int, np.integer],
                       func: typing.Callable[[int, np.ndarray], float], *,
                       funcnr: int = 0, physq: int = 0, degfd: int = 0,
                       f: np.ndarray | None = None) -> np.ndarray | None:
    """
    Fill system vector.

    NOTE: If f is present, this function modifies the input arguments `f`
          in place. If f is not present, f is returned by the function.

    Parameters
    ----------
    mesh : object
        Mesh structure.
    problem : object
        Problem structure.
    geometry : str
        The 'geometry' to fill degrees on:
        - 'nodes': in the nodes given by numbers
        - 'points': in the points given by numbers
        - 'curves': in the curves given by numbers
    numbers : array_like
        An array of 'geometry' numbers.
    func : callable
        Scalar function for filling. Only argument is x, a coordinate vector.

    Keyword arguments
    -----------------
    funcnr : int, optional
        Function number for func, i.e. func(funcnr, x). Default is 1.
    physq : int, optional
        Physical quantity number. Default is 1.
    degfd : int, optional
        Degree of freedom within the physical quantity. Default is 1.
    f : numpy.ndarray, optional
        An existing system vector which will be modified in-place.

    Returns
    -------
    f : numpy.ndarray (optional)
        Filled system vector. Will only be returned if f is not present in the
        function call (if it is present, f is modified in place).

    """

    f_present = f is not None

    if not f_present:
        f = np.zeros(problem.numdegfd)  # type: ignore

    # Convert numbers to a list if an int is supplied
    if isinstance(numbers, (int, np.integer)):
        numbers = typing.cast(list[int], [numbers])

    # Define geometry
    if geometry == 'nodes':
        nodes = numbers
    elif geometry == 'points':
        nodes = mesh.points[numbers]
    elif geometry == 'curves':
        # nnodes = sum([mesh.curves[curve].nnodes for curve in numbers])
        nodes = np.array([], dtype=int)
        for curve in numbers:
            nodes = np.append(nodes, mesh.curves[curve].nodes)
    else:
        raise ValueError(f"Invalid geometry: {geometry}")

    # Fill vector
    for node in nodes:
        posn, ndof = pos_array(problem, node, physq=physq, order='ND')

        if degfd > ndof[0]:
            continue

        f[posn[0][degfd]] = func(funcnr, mesh.coor[node])  # type: ignore

    if not f_present:
        return f  # type: ignore
