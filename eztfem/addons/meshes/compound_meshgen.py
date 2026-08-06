"""Module for compound 2D mesh generation."""
import typing

from ...core.meshgen import quadrilateral2d, Mesh
from ...core.meshgen_extra import mesh_merge

_EltypeLiteral = typing.Literal["tria3", "tria4", "tria6", "tria7",
                                "quad4", "quad9", "quad5"]


def two_blocks2d(
    ne: typing.Sequence[int],
    eltype: _EltypeLiteral,
    *,
    origin: typing.Sequence[float] | None = None,
    length: typing.Sequence[float] | None = None,
    factor: typing.Sequence[float] | None = None,
) -> Mesh:
    """
    Generate mesh for two rectangles side by side.

    Parameters
    ----------
    ne : list of int
        List with three integers [n1, n2, n3], number of elements.
    eltype : str
        Element type (refer to `quadrilateral2d` for details).
    origin : tuple of float, optional
        Coordinates of the "origin" of the domain. Default is [0, 0].
    length : list of float, optional
        Lengths of the domain. Default is [1, 1, 1].
    factor : list of float, optional
        Factor for refinement for each of the three boundaries.
        Default is [1, 1, 1].

    Returns
    -------
    mesh : Mesh
        Mesh object for the domain.

    Notes
    -----
    The domain is divided into two rectangles side by side:

    ::

          L3 -----------------------
             |        |            |
             |        |            |
         n3  |   1    |     2      |
             |        |            |
             |        |            |
           0 -----------------------
            -L1       0            L2
                 n1         n2

    """
    if origin is None:
        origin = (0, 0)
    if length is None:
        length = (1, 1, 1)
    if factor is None:
        factor = (1, 1, 1)

    n1, n2, n3 = ne
    o1, o2 = origin
    l1, l2, l3 = length
    f1, f2, f3 = factor

    mesh1 = quadrilateral2d([n1, n3], eltype, origin=[o1-l1, o2],
                            length=[l1, l3], ratio=[3, 3, 1, 1],
                            factor=[f1, f3, f1, f3])

    mesh2 = quadrilateral2d([n2, n3], eltype, origin=[o1, o2],
                            length=[l2, l3], ratio=[1, 3, 3, 1],
                            factor=[f2, f3, f2, f3])

    mesh = mesh_merge(mesh1, mesh2, curves1=[1], curves2=[3],
                      dir_curves2=[-1], deletecurves1=[1])

    return mesh


def l_shape2d(
    ne: typing.Sequence[int],
    eltype: _EltypeLiteral,
    *,
    origin: typing.Sequence[float] | None = None,
    length: typing.Sequence[float] | None = None,
    factor: typing.Sequence[float] | None = None,
) -> Mesh:
    """
    Generate mesh for an L-shape region (three rectangles).

    Parameters
    ----------
    ne : list of int
        List with four integers [n1, n2, n3, n4], number of elements.
    eltype : str
        Element type (refer to `quadrilateral2d` for details).
    origin : tuple of float, optional
        Coordinates of the "origin" of the domain. Default is (0, 0).
    length : list of float, optional
        Lengths and width of the domain. Default is [1, 1, 1, 1].
    factor : list of float, optional
        Factor for refinement for each of the four boundaries.
        Default is [1, 1, 1, 1].

    Returns
    -------
    mesh : Mesh
        Mesh object for the domain.

    """
    if origin is None:
        origin = (0, 0)
    if length is None:
        length = (1, 1, 1, 1)
    if factor is None:
        factor = (1, 1, 1, 1)

    n1, n2, n3, n4 = ne
    o1, o2 = origin
    l1, l2, l3, l4 = length
    f1, f2, f3, f4 = factor

    mesh1 = quadrilateral2d([n1, n3], eltype, origin=(o1-l1, o2),
                            length=[l1, l3], ratio=[3, 1, 1, 3],
                            factor=[f1, f3, f1, f3])

    mesh2 = quadrilateral2d([n2, n3], eltype, origin=(o1, o2),
                            length=[l2, l3], ratio=[1, 1, 3, 3],
                            factor=[f2, f3, f2, f3])

    mesh3 = mesh_merge(mesh1, mesh2, curves1=[1], curves2=[3],
                       dir_curves2=[-1], deletecurves1=[1])

    mesh4 = quadrilateral2d([n2, n4], eltype, origin=(o1, o2-l4),
                            length=[l2, l4], ratio=[1, 3, 3, 1],
                            factor=[f2, f4, f2, f4])

    mesh = mesh_merge(mesh3, mesh4, curves1=[3], curves2=[2],
                      dir_curves2=[-1], deletecurves1=[3])

    return mesh
