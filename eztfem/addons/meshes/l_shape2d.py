from ...src.quadrilateral2d import quadrilateral2d
from ...src.mesh_merge import mesh_merge

def l_shape2d(ne, eltype, **kwargs):
    """
    Generate mesh for an L-shape region (three rectangles).
    Parameters:
      - ne: list with four integers [n1,n2,n3,n4], number of elements
      - eltype: element type (refer to `quadrilateral2d` for details)
      - origin: tuple with two floats, the coordinates of the "origin" of the domain. Default=(0,0)
      - length: list with four floats, lengths and width of the domain. Default=[1,1,1,1]
      - factor: list with four floats, factor for refinement for each of the four boundaries. Default=[1,1,1,1]
    Returns:
      - mesh: mesh structure for the domain
    """

    # optional arguments
    o = kwargs.get('origin', [0, 0])
    l = kwargs.get('length', [1, 1, 1, 1])
    f = kwargs.get('factor', [1, 1, 1, 1])

    n1, n2, n3, n4 = ne
    o1, o2 = o
    l1, l2, l3, l4 = l
    f1, f2, f3, f4 = f

    mesh1 = quadrilateral2d([n1, n3], eltype, origin=(o1-l1, o2), length=[l1, l3], ratio=[3, 1, 1, 3], factor=[f1, f3, f1, f3])

    mesh2 = quadrilateral2d([n2, n3], eltype, origin=(o1, o2), length=[l2, l3], ratio=[1, 1, 3, 3], factor=[f2, f3, f2, f3])

    mesh3 = mesh_merge(mesh1, mesh2, curves1=[1], curves2=[3], dir_curves2=[-1], deletecurves1=[1])

    mesh4 = quadrilateral2d([n2, n4], eltype, origin=(o1, o2-l4), length=[l2, l4], ratio=[1, 3, 3, 1], factor=[f2, f4, f2, f4])

    mesh = mesh_merge(mesh3, mesh4, curves1=[3], curves2=[2], dir_curves2=[-1], deletecurves1=[3])

    return mesh