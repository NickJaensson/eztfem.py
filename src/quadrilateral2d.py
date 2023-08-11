import numpy as np
from src.mesh_defs import Mesh, Geometry
from src.distribute_elements import distribute_elements

def quadrilateral2d(ne, eltype, **kwargs):
    """"
    Simple meshgenerator for 2D on the region [0,1]x[0,1] with an optional
    translation and/or scaling or deformation of the region.
    input:
      ne=[nx,ny]: number of elements in x and y direction
      eltype: element type
              eltype='tria3': 3 node triangle
              eltype='tria4': 4 node triangle
              eltype='tria6': 6 node triangle
              eltype='tria7': 7 node triangle
              eltype='quad4': 4 node quadrilateral
              eltype='quad9': 9 node quadrilateral
              eltype='quad5': 5 node quadrilateral
    optional arguments:
      string, value couples to set options:
      'origin', origin of the domain [ox,uy]
      'length', length and width of the domain [lx,ly]
      'vertices', give four vertices of the domain in the format
                 [ x1 y1; x2 y2; x3 y3; x4 y4 ]
                note, that origin/length arguments are ignored if vertices are
                given.
      'ratio', a vector [ratio1 ratio2 ratio3 ratio4] where for each of the four
          curves the ratio is given, where
          ratio = 0: equidistant mesh
          ratio = 1: the size of the last element is factor times the first
          ratio = 2: the size of an element is factor times the previous one
          ratio = 3: the size of the last element is 1/factor times the first
          ratio = 4: the size of an element is 1/factor times the previous one
         default=0
      'factor' a vector [factor1,factor2,factor3,factor4] where for each of
          the four curves the factor is given.
         default=[1 1 1 1]
    Examples:
      mesh = quadrilateral2d ( [10,10], 3, 'origin', [1,2], 'length', [2,3] )
    to create a 10x10 mesh of 3 node triangles, where the left lower corner of
    the domain is (1,2) and the size of the domain is 2 by 3.
      mesh = quadrilateral2d ( [10,10], 5, 'vertices', [1,1;3,2;4,5;-1,4] )
    to create a 10x10 mesh of 4 node quadrilateral elements, where the left
    lower corner, right lower corner, upper right corner and upper left corner
    of the domain are (1,1), (3,2), (4,5) and (-1,4), respectively. The edges
    are straight.
    output:
      mesh: mesh structure, having the components
         ndim: dimension of space (ndim=1 or 2)
         nnodes: number of nodes
         coor: array of size (nnodes,ndim), where coor(i,:) are the
               coordinates of node i, with 1<=i<=nnodes.
         nelem: number of elements
         elshape: shape number of the elements. Table of shapes:
               elshape=1: 2 node line elements
               elshape=2: 3 node line elements
               elshape=3: 3 node triangle
               elshape=4: 6 node triangle
               elshape=5: 4 node quadrilateral
               elshape=6: 9 node quadrilateral
               elshape=7: 7 node triangle
               elshape=9: 5 node quadrilateral
               elshape=10: 4 node triangle
         elnumnod: number of nodes in a single element
         topology: array of size (elnumnod,nelem), where topology(:,elem) is
                   an array of global node numbers element elem is connected to
         npoints: number of points
         points: an array containing the node numbers of the points, hence
                 points(i) is the node of point i, with 1<=i<=npoints.
         ncurves: number of curves
         curves: an array of structures, where each structure has the components
             ndim: dimension of space (ndim=2)
             nnodes: number of nodes
             nelem: number of elements
             elshape: shape number of the elements. Table of shapes:
               elshape=1: 2 node elements
               elshape=2: 3 node elements
             elnumnod: number of nodes in a single element
             nodes: array of size nnodes containing the global node numbers,
               i.e. nodes(i) is the global node number of local node i.
             topology: array of size (elnumnod,nelem,2), where
               topology(:,elem,1) is an array of local node numbers element
               elem is connected to.
               topology(:,elem,2) is an array of global node numbers element
               elem is connected to.a
    """

    # optional arguments
    ori = kwargs.get('origin', None)
    length = kwargs.get('length', None)
    def_verts = kwargs.get('vertices', None)
    ratio = kwargs.get('ratio', None)
    factor = kwargs.get('factor', [1.0, 1.0, 1.0, 1.0])
    
    # mesh
    if eltype == 'tria3':
        mesh = rectangle2d_tria3(ne, ratio, factor)
    elif eltype == 'tria4':
        mesh = rectangle2d_tria4(ne, ratio, factor)
    elif eltype == 'tria6':
        mesh = rectangle2d_tria6(ne, ratio, factor)
    elif eltype == 'tria7':
        mesh = rectangle2d_tria7(ne, ratio, factor)
    elif eltype == 'quad4':
        mesh = rectangle2d_quad4(ne, ratio, factor)
    elif eltype == 'quad5':
        mesh = rectangle2d_quad5(ne, ratio, factor)
    elif eltype == 'quad9':
        mesh = rectangle2d_quad9(ne, ratio, factor)
    else:
        raise ValueError(f"Invalid eltype = {eltype}")

    # translate, scale or deform unit region
    if def_verts is not None:

        #     x = [ x1 * (1-xi) + x2 * xi ] (1 - eta) +
        #         [ x4 * (1-xi) + x3 * xi ] eta 
        #     with xi = mesh.coor(1), eta = mesh.coor(2) \in [0,1] from rectangle2d

        coor = np.zeros_like(mesh.coor)
        for i in range(2):
            coor[:, i] = (def_verts[0, i] * (1 - mesh.coor[:, 0]) +
                          def_verts[1, i] * mesh.coor[:, 0]) * (1 - mesh.coor[:, 1]) + \
                          (def_verts[3, i] * (1 - mesh.coor[:, 0]) +
                          def_verts[2, i] * mesh.coor[:, 0]) * mesh.coor[:, 1]
        mesh.coor[:, :] = coor
    else:
        if length is not None:
            mesh.coor *= length
        if ori is not None:
            mesh.coor += ori

    return mesh


def rectangle2d_tria3(ne, ratio, factor):
    # Implementation of rectangle2d_tria3
    raise NotImplementedError


def rectangle2d_tria4(ne, ratio, factor):
    # Implementation of rectangle2d_tria4
    raise NotImplementedError


def rectangle2d_tria6(ne, ratio, factor):
    # Implementation of rectangle2d_tria6
    raise NotImplementedError


def rectangle2d_tria7(ne, ratio, factor):
    # Implementation of rectangle2d_tria7
    raise NotImplementedError


def rectangle2d_quad4(ne, ratio, factor):
    # Implementation of rectangle2d_quad4
    raise NotImplementedError


def rectangle2d_quad5(ne, ratio, factor):
    # Implementation of rectangle2d_quad5
    raise NotImplementedError


def rectangle2d_quad9(ne, ratio, factor):
    """
    Generate a mesh on region [0,1]x[0,1] using quad elements with 9 nodes.
    The numbering for both nodal points and elements is row by row.
    For example a 2x2 mesh is numbered as follows:

          21 -- 22 -- 23 -- 24 -- 25
           |           |           |
          16    17    18    19    20
           |           |           |
          11 -- 12 -- 13 -- 14 -- 15
           |           |           |
           6     7     8     9    10
           |           |           |
           1 --  2 --  3 --  4 --  5

       Elements:

           x --  x --  x --  x --  x
           |           |           |
           x     3     x     4     x
           |           |           |
           x --  x --  x --  x --  x
           |           |           |
           x     1     x     2     x
           |           |           |
           x --  x --  x --  x --  x

      Also generated are four points and four curves:

                   <----- 
                      C3
          P4 --------------------- P3
           |                       |           Note:   C1=P1-P2
           |                       |                   C2=P2-P3
      |    |                       |    ^              C3=P3-P4
      | C4 |                       | C2 |              C4=P4-P1
     \|/   |                       |    |
           |                       |
           |                       |
          P1 --------------------- P2
                      C1
                    ----->
    """
    
    print('rectangle_quad9')

    n = ne[0]
    m = ne[1]

    numnod = 9 # number of nodes per element
    nn1row = 2 * n + 1 # number of nodes in one row
    nn1col = 2 * m + 1 # number of nodes in one column

    # create mesh and fill structure
    mesh = Mesh()
    mesh.ndim = 2
    mesh.nnodes = nn1row * nn1col
    mesh.elshape = 6
    mesh.nelem = n * m
    mesh.elnumnod = numnod
    mesh.npoints = 4
    mesh.ncurves = 4

    mesh.topology = np.zeros((numnod, mesh.nelem),dtype=int)
    mesh.coor = np.zeros((mesh.nnodes, mesh.ndim))
    mesh.points = np.zeros((mesh.npoints),dtype=int)

    # topology
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            nn1 = 2 * (j - 1) * nn1row + 2 * (i - 1) # number of nodes before element (i,j)
            nn2 = nn1 + nn1row # number of nodes before element (i,j) + 1 row
            nn3 = nn2 + nn1row # number of nodes before element (i,j) + 2 rows

            e = i + (j - 1) * n

            mesh.topology[0, e - 1] = nn1 + 1
            mesh.topology[1, e - 1] = nn1 + 2
            mesh.topology[2, e - 1] = nn1 + 3
            mesh.topology[3, e - 1] = nn2 + 3
            mesh.topology[4, e - 1] = nn3 + 3
            mesh.topology[5, e - 1] = nn3 + 2
            mesh.topology[6, e - 1] = nn3 + 1
            mesh.topology[7, e - 1] = nn2 + 1
            mesh.topology[8, e - 1] = nn2 + 2

    mesh.topology = mesh.topology - 1

    # coordinates
    if ratio is None:
        # equidistant
        deltax = 0.5 / n
        deltay = 0.5 / m
        for i in range(1, nn1row + 1):
            for j in range(1, nn1col + 1):
                node = i + (j - 1) * nn1row - 1
                mesh.coor[node, 0] = (i - 1) * deltax
                mesh.coor[node, 1] = (j - 1) * deltay
    else:
        # non-equidistant
        x1 = np.zeros(nn1row)
        x2 = np.zeros(nn1col)
        x3 = np.zeros(nn1row)
        x4 = np.zeros(nn1col)
        x1[0::2] = distribute_elements(n, ratio[0], factor[0])
        x2[0::2] = distribute_elements(m, ratio[1], factor[1])
        x3[0::2] = distribute_elements(n, ratio[2], factor[2])
        x4[0::2] = distribute_elements(m, ratio[3], factor[3])

        # mid-side nodes
        for elem in range(1, n + 1):
            k = 2 * elem - 2
            x1[k + 1] = (x1[k] + x1[k + 2]) / 2
            x3[k + 1] = (x3[k] + x3[k + 2]) / 2

        for elem in range(1, m + 1):
            k = 2 * elem - 2
            x2[k + 1] = (x2[k] + x2[k + 2]) / 2
            x4[k + 1] = (x4[k] + x4[k + 2]) / 2

        x3 = 1 - x3[::-1]
        x4 = 1 - x4[::-1]

        # create straight lines in reference square [0,1]x[0,1]
        for i in range(1, nn1row + 1):
            for j in range(1, nn1col + 1):
                node = i + (j - 1) * nn1row - 1
                D = (x1[i - 1] - x3[i - 1]) * (x4[j - 1] - x2[j - 1]) - 1
                mesh.coor[node, 0] = (x4[j - 1] * (x1[i - 1] - x3[i - 1]) - x1[i - 1]) / D
                mesh.coor[node, 1] = (x1[i - 1] * (x4[j - 1] - x2[j - 1]) - x4[j - 1]) / D

    # points
    mesh.points[0] = 0
    mesh.points[1] = nn1row - 1
    mesh.points[2] = nn1row * nn1col - 1
    mesh.points[3] = nn1row * (nn1col - 1)

    # curves
    a = [nn1row, nn1col, nn1row, nn1col]
    b = [n, m, n, m]

    mesh.curves = [Geometry(elshape=2,ndim=2,elnumnod=3,nnodes=a[curve],nelem=b[curve]) for curve in range(4)]

    mesh.curves[0].nodes = np.arange(1, nn1row + 1) - 1
    mesh.curves[1].nodes = nn1row * np.arange(1, nn1col + 1, 1) - 1
    mesh.curves[2].nodes = (nn1col - 1) * nn1row + np.arange(nn1row, 0, -1) - 1
    mesh.curves[3].nodes = 1 - nn1row + nn1row * np.arange(nn1col, 0, -1) - 1

    for curve in range(mesh.ncurves):

        nelem = mesh.curves[curve].nelem
        elnumnod = mesh.curves[curve].elnumnod
        mesh.curves[curve].topology = np.zeros((elnumnod, nelem, 2), dtype=int)

        # local numbering
        for elem in range(nelem):
            mesh.curves[curve].topology[:, elem, 0] = np.arange(2 * elem, 2 * elem + 3)

        # global numbering
        mesh.curves[curve].topology[:, :, 1] = mesh.curves[curve].nodes[ mesh.curves[curve].topology[:, :, 0] ]

    return mesh