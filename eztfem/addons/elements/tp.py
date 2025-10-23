import typing

import numpy as np

from ...core.tp import UserProtocol

P = typing.ParamSpec("P")
T = typing.TypeVar("T")


class CoorsysAwareUser(UserProtocol, typing.Protocol):
    """User protocol that is aware of its coordinate system."""

    coorsys: int


class FuncAwareUser(UserProtocol, typing.Protocol[P, T]):
    """User protocol that can provide user functions to element routines."""

    funcnr: int
    func: typing.Callable[P, T]


# XXX: streamfunction_elements.py expects user.v
class SFVelocityAwareUser(UserProtocol, typing.Protocol):
    """User protocol that is aware of its stream velocity."""

    v: np.ndarray


# XXX: stokes_elements.py and poisson_elements.py expect user.u
class VelocityAwareUser(UserProtocol, typing.Protocol):
    """User protocol that is aware of its stream velocity."""

    u: np.ndarray


class DynamicViscosityAwareUser(UserProtocol, typing.Protocol):
    """User protocol that is aware of its stream dynamic viscosity."""

    mu: float


class AlphaAwareUser(UserProtocol, typing.Protocol):
    """User protocol that is aware of a Poisson flow alpha constant."""

    alpha: float


class PoissonElemUser(AlphaAwareUser, CoorsysAwareUser,
                      FuncAwareUser[[int, np.ndarray], np.ndarray]):
    """User protocol that provides the fields required for
    :func:`poisson_elem`."""
    pass


class StokesElemUser(CoorsysAwareUser, DynamicViscosityAwareUser,
                     FuncAwareUser[[int, np.ndarray], np.ndarray],
                     VelocityAwareUser, typing.Protocol):
    """User protocol that provides the fields required for
    :func:`stokes_elem`."""
    pass


class StokesDerivUser(CoorsysAwareUser, VelocityAwareUser,
                      typing.Protocol):
    """User protocol that provides the fields required for
    :func:`stokes_deriv`."""

    # TODO: Maybe these can be clarified more.
    comp: typing.Literal[0, 1, 2, 3, 4, 5, 6, 7]
    """The component to derive:
    
    0. dudx,
    1. dudy,
    2. dvdx,
    3. dvdy,
    4. Vorticity,
    5. Gradient of velocity field in theta direction,
    6. Divergence of velocity field,
    7. Effective strain rate.
    """


class NatbounCurveUser(CoorsysAwareUser,
                             FuncAwareUser[[int, np.ndarray], np.ndarray],
                             typing.Protocol):
    """User protocol that provides the fields required for
    :func:`poisson_natboun_curve` and :func:`stokes_natboun_curve`."""
    pass


class StokesFlowrateCurveUser(CoorsysAwareUser, VelocityAwareUser,
                              typing.Protocol):
    """User protocol that provides the fields required for
    :func:`stokes_flowrate_curve`."""
    pass


class StreamfunctionUser(CoorsysAwareUser, SFVelocityAwareUser,
                          typing.Protocol):
    """User protocol that provides the fields required for
    :func:`streamfunction_elem` and :func:`streamfunction_natboun_curve`."""
    pass

