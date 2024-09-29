from ...core.quadrilateral2d import quadrilateral2d
from ...core.mesh_merge import mesh_merge


def two_blocks2d(ne, eltype, **kwargs):
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

    # Optional arguments
    oo = kwargs.get('origin', [0, 0])
    ll = kwargs.get('length', [1, 1, 1])
    ff = kwargs.get('factor', [1, 1, 1])

    n1, n2, n3 = ne
    o1, o2 = oo
    l1, l2, l3 = ll
    f1, f2, f3 = ff

    mesh1 = quadrilateral2d([n1, n3], eltype, origin=[o1-l1, o2],
                            length=[l1, l3], ratio=[3, 3, 1, 1],
                            factor=[f1, f3, f1, f3])

    mesh2 = quadrilateral2d([n2, n3], eltype, origin=[o1, o2],
                            length=[l2, l3], ratio=[1, 3, 3, 1],
                            factor=[f2, f3, f2, f3])

    mesh = mesh_merge(mesh1, mesh2, curves1=[1], curves2=[3],
                      dir_curves2=[-1], deletecurves1=[1])

    return mesh
