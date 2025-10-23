import typing

import numpy as np
import typing_extensions


py_dtypeT = typing_extensions.TypeVar("py_dtypeT", default=typing.Any)
np_dtypeT = typing_extensions.TypeVar("np_dtypeT", bound=np.generic, default=typing.Any)

Array1D: typing.TypeAlias = np.ndarray[tuple[int], np.dtype[np_dtypeT]]
Array2D: typing.TypeAlias = np.ndarray[tuple[int, int], np.dtype[np_dtypeT]]
ArrayLike: typing.TypeAlias = (
    py_dtypeT  # A python type, e.g. int
    | np_dtypeT  # An equivalent numpy type, e.g. np.integer
    | typing.Sequence[py_dtypeT]  # Any ordered sequence of the python type, e.g. list[int]
    | typing.Sequence[np_dtypeT]  # Any ordered sequence of the numpy type, e.g. list[np.integer]
    | np.typing.NDArray[np_dtypeT]  # A numpy array of the numpy type and arbitrary shape
)

IntArrayLike: typing.TypeAlias = ArrayLike[int, np.integer]

Order: typing.TypeAlias = typing.Literal["DN", "ND"]
Ratio: typing.TypeAlias = typing.Literal[0, 1, 2, 3, 4]

NestedSequence: typing.TypeAlias = typing.Sequence[typing.Sequence[py_dtypeT]]


class User(typing.Protocol):
    """Protocol class describing a minimal User class."""

    wg: np.ndarray
    xr: np.ndarray
    phi: np.ndarray
    dphi: np.ndarray
    psi: np.ndarray


RoutineReturnT = typing.TypeVar("RoutineReturnT")

ElementRoutine: typing.TypeAlias = typing.Callable[  # A function that takes
    [int, np.ndarray, User, NestedSequence[int]],  # elem, coor, user, pos 
    RoutineReturnT,  # and returns a yet-to-be-specified type
]

PosvecElementRoutine: typing.TypeAlias = typing.Callable[
    [int, np.ndarray, User, NestedSequence[int], NestedSequence[int]],  # elem, coor, user, pos, posvec
    RoutineReturnT,
]
