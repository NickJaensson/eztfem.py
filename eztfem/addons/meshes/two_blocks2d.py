from eztfem.src.quadrilateral2d import quadrilateral2d
from eztfem.src.mesh_merge import mesh_merge

def two_blocks2d(ne, eltype, **kwargs):
    """
    TWO_BLOCKS2D  generate mesh for two rectangles side by side
      mesh = TWO_BLOCKS2D ( ne, eltype, 'option1', value1, .... )
      input:
        ne=[n1,n2,n3]: number of elements (see below)
        eltype: element type (see quadrilateral2d)
      optional arguments:
        string, value couples to set options:
        'origin', the coordinates of the "origin" of the domain [ox,oy] 
           default=[0 0]
        'length', lengths of the domain [L1,L2,L3]
           default=[1 1 1]
        'factor' a vector [factor1,factor2,factor3] where for each of the three
                 boundaries the factor for refinement is given.
           default=[1 1 1]
      output:
        mesh: mesh structure for the domain:
    
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
    o = [0, 0]
    l = [1, 1, 1]
    f = [1, 1, 1]

    if 'origin' in kwargs:
        o = kwargs['origin']
    if 'length' in kwargs:
        l = kwargs['length']
    if 'factor' in kwargs:
        f = kwargs['factor']

    n1, n2, n3 = ne
    o1, o2 = o
    l1, l2, l3 = l
    f1, f2, f3 = f

    mesh1 = quadrilateral2d ( [n1, n3], eltype, origin=[o1-l1,o2], length=[l1,l3], ratio=[3,3,1,1], factor=[f1,f3,f1,f3] )
  
    mesh2 = quadrilateral2d ( [n2, n3], eltype, origin=[o1,o2],length=[l2,l3], ratio=[1,3,3,1], factor=[f2,f3,f2,f3] )
  
    mesh = mesh_merge ( mesh1, mesh2, curves1=[1], curves2=[-3], deletecurves1=[1] )

    return mesh