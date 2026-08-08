"""Module to define the User class"""
import typing

import numpy as np
import numpy.typing as npt

FloatArray: typing.TypeAlias = npt.NDArray[np.floating]

# Element routines (in eztfem.addons.elements and user code) all share the
# call signature (elem, coor, user, pos[, posvec]), but their return shape
# varies by call site (build_system expects (elemmat, elemvec),
# deriv_vector expects elemvec alone, etc.), so this is intentionally
# left as a loose Callable rather than a single precise Protocol.
ElementRoutine: typing.TypeAlias = typing.Callable[..., typing.Any]

# User.func is called as func(funcnr, x) by element/derivative routines, but
# its argument count, and whether it returns a scalar or a vector, varies by
# problem (see e.g. examples/poisson/func.py vs examples/stokes/traction_func
# .py), so this is intentionally left as a loose Callable, matching
# ElementRoutine above.
UserFunc: typing.TypeAlias = typing.Callable[..., typing.Any]


class User:
    """
    Define the User class: this class is mainly a collection of parameters and
    data that should be accessible to the element routines.

    Attributes
    ----------
    wg : np.ndarray
        Weights of Gauss points.
    xr : np.ndarray
        Location of Gauss points.
    phi : np.ndarray
        Values of phi shape function in Gauss points.
    dphi : np.ndarray
        Values of phi shape function derivatives in Gauss points.
    psi : np.ndarray
        Values of psi shape function in Gauss points.
    coorsys : int
        Coordinate system selector used by element routines (set by user
        code, not by __init__).
    alpha : float
        Diffusion coefficient used by the Poisson element routines (set by
        user code, not by __init__).
    mu : float
        Viscosity used by the Stokes element routines (set by user code,
        not by __init__).
    funcnr : int
        Selector passed as the first argument to `func` (set by user code,
        not by __init__).
    func : callable
        Right-hand-side/boundary function, called as func(funcnr, x) (set by
        user code, not by __init__).
    comp : int
        Selects which derived quantity `stokes_deriv` computes (set by user
        code, not by __init__).
    u : np.ndarray
        Solution vector used by the `*_deriv` routines (set by user code,
        not by __init__).
    v : np.ndarray
        Velocity vector used by the streamfunction routines (set by user
        code, not by __init__).

    """

    # The attributes below are set by user code after construction (their
    # presence/absence is what User.__eq__ compares), so they are declared
    # here without a value rather than assigned in __init__.
    coorsys: int
    alpha: float
    mu: float
    funcnr: int
    func: UserFunc
    comp: int
    u: FloatArray
    v: FloatArray

    def __init__(self) -> None:
        """
        Initializes the Problem object with the all attributes empty numpy
        arrays (np.array([])).

        """
        self.wg = np.array([])
        self.xr = np.array([])
        self.phi = np.array([])
        self.dphi = np.array([])
        self.psi = np.array([])

    def __eq__(self, other: object) -> bool:
        """
        Checks equivalence of two User objects (overloads == sign).

        Parameters
        ----------
        other : User
            The other User object to compare with.

        Returns
        -------
        bool
            True if the User objects are equivalent, False otherwise.

        Notes
        -----
        See NOTE_ON_COMPARING_ARRAYS.md for the use of np.squeeze.
        Only a selection of possible attributes is checked.

        """
        if not isinstance(other, User):
            return NotImplemented

        # Get list of names of attributes
        attributes_self = list(self.__dict__)
        attributes_other = list(other.__dict__)

        # Check if attributes exist in both Users
        check0 = [
            ('wg' in attributes_self) is ('wg' in attributes_other),
            ('xr' in attributes_self) is ('xr' in attributes_other),
            ('phi' in attributes_self) is ('phi' in attributes_other),
            ('psi' in attributes_self) is ('psi' in attributes_other),
            ('dphi' in attributes_self) is ('dphi' in attributes_other),
            ('coorsys' in attributes_self) is ('coorsys' in attributes_other),
            ('alpha' in attributes_self) is ('alpha' in attributes_other),
            ('mu' in attributes_self) is ('mu' in attributes_other),
            ('funcnr' in attributes_self) is ('funcnr' in attributes_other)
        ]

        if not all(check0):
            print("WARNING: Users do not have equal attributes:")
            print(check0)
            return False

        # Check if existing attributes have same values
        check1 = [True]  # Avoid empty list
        if 'wg' in attributes_self:
            check1.append(np.allclose(self.wg, other.wg, atol=1e-15, rtol=0))
        if 'xr' in attributes_self:
            check1.append(np.allclose(np.squeeze(self.xr), other.xr,
                                      atol=1e-15, rtol=0))
        if 'phi' in attributes_self:
            check1.append(np.allclose(np.squeeze(self.phi), other.phi,
                                      atol=1e-15, rtol=0))
        if 'psi' in attributes_self:
            check1.append(np.allclose(np.squeeze(self.psi), other.psi,
                                      atol=1e-15, rtol=0))
        if 'dphi' in attributes_self:
            check1.append(np.allclose(self.dphi, other.dphi, atol=1e-15,
                                      rtol=0))
        if 'coorsys' in attributes_self:
            check1.append(self.coorsys == other.coorsys)
        if 'alpha' in attributes_self:
            check1.append(self.alpha == other.alpha)
        if 'mu' in attributes_self:
            check1.append(self.mu == other.mu)
        if 'funcnr' in attributes_self:
            check1.append(self.funcnr == other.funcnr)

        # Print a warning when not equivalent (for debugging purposes)
        if not all(check1):
            print("WARNING: User attributes not equivalent:")
            print(check1)

        return all(check1)
