import numpy as np
from src.mesh_class import Mesh, Geometry
from src.distribute_elements import distribute_elements

def quadrilateral2d(num_el, eltype, **kwargs):
    """"
    QUADRILATERAL2D  Simple meshgenerator for quadrilateral 2D regions.
      Simple meshgenerator for 2D on the region [0,1]x[0,1] with an optional
      translation and/or scaling or deformation of the region.
      mesh = QUADRILATERAL2D ( ne, eltype, 'option1', value1, .... )
      input:
        num_el=[nx,ny]: number of elements in x and y direction
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
        'ratio', a vector [ratio1 ratio2 ratio3 ratio4] where for each of the 
            four curves the ratio is given, where
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
                     an array of global node numbers element elem is connected 
                     to
           npoints: number of points
           points: an array containing the node numbers of the points, hence
                   points(i) is the node of point i, with 1<=i<=npoints.
           ncurves: number of curves
           curves: an array of structures, where each structure has the 
               components
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
                 elem is connected to.
             
    """

    # optional arguments
    ori = kwargs.get('origin', None)
    length = kwargs.get('length', None)
    def_verts = kwargs.get('vertices', None)
    ratio = kwargs.get('ratio', None)
    factor = kwargs.get('factor', None)

    # mesh
    if eltype == 'tria3': # 3 node triangle
        mesh = rectangle2d_tria3(num_el, ratio, factor)
    elif eltype == 'tria4': # 4 node triangle
        mesh = rectangle2d_tria4(num_el, ratio, factor)
    elif eltype == 'tria6': # 6 node triangle
        mesh = rectangle2d_tria6(num_el, ratio, factor)
    elif eltype == 'tria7': # 7 node triangle
        mesh = rectangle2d_tria7(num_el, ratio, factor)
    elif eltype == 'quad4': # 4 node quadrilateral
        mesh = rectangle2d_quad4(num_el, ratio, factor)
    elif eltype == 'quad5': # 5 node quadrilateral
        mesh = rectangle2d_quad5(num_el, ratio, factor)
    elif eltype == 'quad9': # 9 node quadrilateral
        mesh = rectangle2d_quad9(num_el, ratio, factor)
    else:
        raise ValueError(f"Invalid eltype = {eltype}")

    # translate, scale or deform unit region
    if def_verts is not None:
        # x = [ x1 * (1-xi) + x2 * xi ] (1 - eta) +
        #     [ x4 * (1-xi) + x3 * xi ] eta
        # with xi = mesh.coor(1), eta = mesh.coor(2) \in [0,1] from rectangle2d
        x_coor = mesh.coor[:, 0]
        y_coor = mesh.coor[:, 1]
        x_weights_1 = (1 - x_coor)[:, np.newaxis]
        x_weights_2 = x_coor[:, np.newaxis]
        y_weights_1 = (1 - y_coor)[:, np.newaxis]
        y_weights_2 = y_coor[:, np.newaxis]
        coor_x_1 = x_weights_1 * def_verts[0, :] + x_weights_2 * def_verts[1, :]
        coor_x_2 = x_weights_1 * def_verts[3, :] + x_weights_2 * def_verts[2, :]
        mesh.coor = y_weights_1 * coor_x_1 + y_weights_2 * coor_x_2
    else:
        if length is not None:
            mesh.coor *= length
        if ori is not None:
            mesh.coor += ori

    return mesh

def rectangle2d_tria3(num_el, ratio, factor):
    """
    Generate a mesh on region [0,1]x[0,1] using triangular elements with 3 nodes.
    The numbering for both nodal points and elements is row by row.
    For example a 2x2 mesh is numbered as follows
   
      Nodes:
   
          7 ------- 8 ------- 9
          |       / |       / |
          |    /    |    /    |
          | /       | /       |
          4 ------- 5 ------- 6
          |       / |       / |
          |    /    |    /    |
          | /       | /       |
          1 ------- 2 ------- 3
   
      Elements:
   
          x ------- x ------- x
          |  5    / |  7    / |
          |    /    |    /    |
          | /     6 | /    8  |
          x ------- x ------- x
          |  1    / |  3    / |
          |    /    |    /    |
          | /    2  | /    4  |
          x ------- x ------- x
   
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
    
    print('rectangle2d_tria3')

    n_x, n_y = num_el

    numnod = 3 # number of nodes per element
    nn1row = n_x + 1 # number of nodes in one row
    nn1col = n_y + 1 # number of nodes in one column

    # create mesh and fill structure
    nelem = 2 * n_x * n_y
    nnodes = nn1row * nn1col
    mesh = Mesh(ndim=2, nnodes=nnodes, elshape=3, nelem=nelem, 
                elnumnod=numnod, npoints=4, ncurves=4, 
                topology=np.zeros((numnod, nelem),dtype=int),
                coor=np.zeros((nnodes, 2)), points=np.zeros((4),dtype=int))

    # topology
    for i in range(n_x):
        for j in range(n_y):
            nn1 = j * nn1row + i # number of nodes before element (i,j)
            nn2 = nn1 + nn1row # number of nodes before element (i,j) + 1 row

            elem = 2 * i + 2 * j * n_x # element number
            mesh.topology[:, elem] = [nn1, nn2+1, nn2]
            elem += 1 # element number
            mesh.topology[:, elem] = [nn1, nn1+1, nn2+1]

    # coordinates
    if ratio is None:
        # equidistant
        deltax = 1.0 / n_x
        deltay = 1.0 / n_y
        for i in range(nn1row):
            for j in range(nn1col):
                node = i + j * nn1row
                mesh.coor[node] = [i * deltax, j * deltay]
    else:
        # non-equidistant
        x1, x3 = np.zeros(nn1row), np.zeros(nn1row),
        x2, x4 = np.zeros(nn1col), np.zeros(nn1col)
        arrays = [x1, x2, x3, x4]
        for idx, array in enumerate(arrays):
            array[:] = distribute_elements(n_x if idx % 2 == 0 
                                             else n_y, ratio[idx], factor[idx])

        x3 = 1 - x3[::-1]
        x4 = 1 - x4[::-1]

        # create straight lines in reference square [0,1]x[0,1]
        i_vals = np.arange(nn1row)
        j_vals = np.arange(nn1col)

        I, J = np.meshgrid(i_vals, j_vals)

        nodes = I + J * nn1row
        D = (x1[I] - x3[I]) * (x4[J] - x2[J]) - 1

        mesh.coor[nodes, 0] = (x4[J] * (x1[I] - x3[I]) - x1[I]) / D
        mesh.coor[nodes, 1] = (x1[I] * (x4[J] - x2[J]) - x4[J]) / D

    # points
    mesh.points = np.array([0, nn1row-1, 
                            nn1row*nn1col-1, nn1row*(nn1col-1)],dtype=int)

    # curves
    temp_a = [nn1row, nn1col, nn1row, nn1col]
    temp_b = [n_x, n_y, n_x, n_y]

    mesh.curves = [Geometry(elshape=1,ndim=2,elnumnod=2,nnodes=temp_a[curve],
                            nelem=temp_b[curve]) for curve in range(4)]

    # nodes of curves
    mesh.curves[0].nodes = np.arange(nn1row)
    mesh.curves[1].nodes = nn1row * ( np.arange(nn1col) + 1 ) - 1
    mesh.curves[2].nodes = (nn1col - 1) * nn1row + np.arange(nn1row)[::-1]
    mesh.curves[3].nodes = nn1row * np.arange(nn1col)[::-1]

    # topology of elements on curves
    for curve in mesh.curves:

        nelem = curve.nelem
        elnumnod = curve.elnumnod
        curve.topology = np.zeros((elnumnod, nelem, 2), dtype=int)

        # local numbering
        for elem in range(nelem):
            curve.topology[:, elem, 0] = np.arange(elem, elem + 2)

        # global numbering
        curve.topology[:, :, 1] = curve.nodes[ curve.topology[:, :, 0] ]

    return mesh

def rectangle2d_tria4(ne, ratio, factor):
    raise NotImplementedError

def rectangle2d_tria6(ne, ratio, factor):
    raise NotImplementedError

def rectangle2d_tria7(ne, ratio, factor):
    raise NotImplementedError

def rectangle2d_quad4(num_el, ratio, factor):
    """
    Generate a mesh on region [0,1]x[0,1] using quad elements with 4 nodes.
    The numbering for both nodal points and elements is row by row.
    For example a 2x2 mesh is numbered as follows    
      Nodes:    
          7 ------- 8 ------- 9
          |         |         |
          |         |         |
          |         |         |
          4 ------- 5 ------- 6
          |         |         |
          |         |         |
          |         |         |
          1 ------- 2 ------- 3    
      Elements:    
          x ------- x ------- x
          |         |         |
          |    3    |    4    |
          |         |         |
          x ------- x ------- x
          |         |         |
          |    1    |    2    |
          |         |         |
          x ------- x ------- x    
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
    
    print('rectangle_quad4')

    n_x, n_y = num_el

    numnod = 4 # number of nodes per element
    nn1row = n_x + 1 # number of nodes in one row
    nn1col = n_y + 1 # number of nodes in one column

    # create mesh and fill structure
    nelem = n_x * n_y
    nnodes = nn1row * nn1col
    mesh = Mesh(ndim=2, nnodes=nnodes, elshape=5, nelem=nelem, 
                elnumnod=numnod, npoints=4, ncurves=4, 
                topology=np.zeros((numnod, nelem),dtype=int),
                coor=np.zeros((nnodes, 2)), points=np.zeros((4),dtype=int))

    # topology
    for i in range(n_x):
        for j in range(n_y):
            nn1 = j * nn1row + i # number of nodes before element (i,j)
            nn2 = nn1 + nn1row # number of nodes before element (i,j) + 1 row

            elem = i + j * n_x # element number

            mesh.topology[:, elem] = [nn1, nn1+1, nn2+1, nn2]

    # coordinates
    if ratio is None:
        # equidistant
        deltax = 1.0 / n_x
        deltay = 1.0 / n_y
        for i in range(nn1row):
            for j in range(nn1col):
                node = i + j * nn1row
                mesh.coor[node] = [i * deltax, j * deltay]
    else:
        # non-equidistant
        x1, x3 = np.zeros(nn1row), np.zeros(nn1row),
        x2, x4 = np.zeros(nn1col), np.zeros(nn1col)
        arrays = [x1, x2, x3, x4]
        for idx, array in enumerate(arrays):
            array[:] = distribute_elements(n_x if idx % 2 == 0 
                                             else n_y, ratio[idx], factor[idx])

        x3 = 1 - x3[::-1]
        x4 = 1 - x4[::-1]

        # create straight lines in reference square [0,1]x[0,1]
        i_vals = np.arange(nn1row)
        j_vals = np.arange(nn1col)

        I, J = np.meshgrid(i_vals, j_vals)

        nodes = I + J * nn1row
        D = (x1[I] - x3[I]) * (x4[J] - x2[J]) - 1

        mesh.coor[nodes, 0] = (x4[J] * (x1[I] - x3[I]) - x1[I]) / D
        mesh.coor[nodes, 1] = (x1[I] * (x4[J] - x2[J]) - x4[J]) / D

    # points
    mesh.points = np.array([0, nn1row-1, 
                            nn1row*nn1col-1, nn1row*(nn1col-1)],dtype=int)

    # curves
    temp_a = [nn1row, nn1col, nn1row, nn1col]
    temp_b = [n_x, n_y, n_x, n_y]

    mesh.curves = [Geometry(elshape=1,ndim=2,elnumnod=2,nnodes=temp_a[curve],
                            nelem=temp_b[curve]) for curve in range(4)]

    # nodes of curves
    mesh.curves[0].nodes = np.arange(nn1row)
    mesh.curves[1].nodes = nn1row * ( np.arange(nn1col) + 1 ) - 1
    mesh.curves[2].nodes = (nn1col - 1) * nn1row + np.arange(nn1row)[::-1]
    mesh.curves[3].nodes = nn1row * np.arange(nn1col)[::-1]

    # topology of elements on curves
    for curve in mesh.curves:

        nelem = curve.nelem
        elnumnod = curve.elnumnod
        curve.topology = np.zeros((elnumnod, nelem, 2), dtype=int)

        # local numbering
        for elem in range(nelem):
            curve.topology[:, elem, 0] = np.arange(elem, elem + 2)

        # global numbering
        curve.topology[:, :, 1] = curve.nodes[ curve.topology[:, :, 0] ]

    return mesh


def rectangle2d_quad5(num_el, ratio, factor):
    """
    Generate a mesh on region [0,1]x[0,1] using quad elements with 5 nodes.
    The numbering for both nodal points and elements is row by row.
    For example a 2x2 mesh is numbered as follows    
      Nodes:    
          11------- 12------- 13
          |         |         |
          |    9    |   10    |
          |         |         |
          6 ------- 7 ------- 8
          |         |         |
          |    4    |    5    |
          |         |         |
          1 ------- 2 ------- 3    
      Elements:    
          x ------- x ------- x
          |         |         |
          |    3    |    4    |
          |         |         |
          x ------- x ------- x
          |         |         |
          |    1    |    2    |
          |         |         |
          x ------- x ------- x    
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
    
    print('rectangle_quad5')

    n_x, n_y = num_el

    numnod = 5 # number of nodes per element
    nn1row = n_x + 1 # number of nodes in one row
    nn1col = n_y + 1 # number of nodes in one column

    # create mesh and fill structure
    nelem = n_x * n_y
    nnodes = nn1row * nn1col + n_x * n_y
    mesh = Mesh(ndim=2, nnodes=nnodes, elshape=9, nelem=nelem, 
                elnumnod=numnod, npoints=4, ncurves=4, 
                topology=np.zeros((numnod, nelem),dtype=int),
                coor=np.zeros((nnodes, 2)), points=np.zeros((4),dtype=int))

    # topology
    for i in range(n_x):
        for j in range(n_y):
            nn1 = j * ( nn1row + n_x ) + i # number of nodes before element (i,j)
            nn2 = nn1 + nn1row # number of nodes before element (i,j) + 1 row
            nn3 = nn2 + n_x # number of nodes before element (i,j) + 2 rows

            elem = i + j * n_x # element number

            mesh.topology[:, elem] = [nn1, nn1+1, nn3+1, nn3, nn2]

    # coordinates
    if ratio is None:
        # equidistant
        deltax = 1.0 / n_x
        deltay = 1.0 / n_y
        for i in range(nn1row):
            for j in range(nn1col):
                node = i + j * ( nn1row + n_x )
                mesh.coor[node] = [i * deltax, j * deltay]
        for i in range(n_x):
            for j in range(n_y):
                node = i + j * ( nn1row + n_x ) + nn1row
                mesh.coor[node] = [i * deltax + deltax/2, j * deltay + deltay/2]
    else:
        # non-equidistant
        x1, x3 = np.zeros(nn1row), np.zeros(nn1row),
        x2, x4 = np.zeros(nn1col), np.zeros(nn1col)
        arrays = [x1, x2, x3, x4]
        for idx, array in enumerate(arrays):
            array[:] = distribute_elements(n_x if idx % 2 == 0 
                                             else n_y, ratio[idx], factor[idx])

        x3 = 1 - x3[::-1]
        x4 = 1 - x4[::-1]

        # create straight lines in reference square [0,1]x[0,1]
        i_vals = np.arange(nn1row)
        j_vals = np.arange(nn1col)

        I, J = np.meshgrid(i_vals, j_vals)

        nodes = I + J *  ( nn1row + n_x )
        D = (x1[I] - x3[I]) * (x4[J] - x2[J]) - 1

        mesh.coor[nodes, 0] = (x4[J] * (x1[I] - x3[I]) - x1[I]) / D
        mesh.coor[nodes, 1] = (x1[I] * (x4[J] - x2[J]) - x4[J]) / D

        i_vals_inner = np.arange(n_x)
        j_vals_inner = np.arange(n_y)

        I_inner, J_inner = np.meshgrid(i_vals_inner, j_vals_inner, 
                                       indexing='ij')

        node = I_inner + J_inner * (nn1row + n_x) + nn1row
        node1 = I_inner + J_inner * (nn1row + n_x)
        node4 = I_inner + (J_inner + 1) * (nn1row + n_x)

        mesh.coor[node,:] = (mesh.coor[node1,:] + mesh.coor[node1+1,:] +
                             mesh.coor[node4,:] + mesh.coor[node4+1,:]) / 4

    # points
    mesh.points = np.array([0, nn1row-1, 
                  nn1row*nn1col+n_x*n_y-1, nn1row*(nn1col-1)+n_x*n_y],dtype=int)

    # curves
    temp_a = [nn1row, nn1col, nn1row, nn1col]
    temp_b = [n_x, n_y, n_x, n_y]

    mesh.curves = [Geometry(elshape=1,ndim=2,elnumnod=2,nnodes=temp_a[curve],
                            nelem=temp_b[curve]) for curve in range(4)]

   # nodes of curves
    mesh.curves[0].nodes = np.arange(1,nn1row+1) - 1
    mesh.curves[1].nodes = ( nn1row + n_x ) * ( np.arange(1,nn1col+1) ) - n_x - 1
    mesh.curves[2].nodes = (nn1col - 1) * nn1row  + n_x * n_y + np.arange(1,nn1row+1)[::-1] - 1
    mesh.curves[3].nodes = 1 - nn1row - n_x + ( nn1row + n_x ) * np.arange(1,nn1col+1)[::-1] - 1

    # topology of elements on curves
    for curve in mesh.curves:

        nelem = curve.nelem
        elnumnod = curve.elnumnod
        curve.topology = np.zeros((elnumnod, nelem, 2), dtype=int)

        # local numbering
        for elem in range(nelem):
            curve.topology[:, elem, 0] = np.arange(elem, elem + 2)

        # global numbering
        curve.topology[:, :, 1] = curve.nodes[ curve.topology[:, :, 0] ]

    return mesh


def rectangle2d_quad9(num_el, ratio, factor):
    """
    Generate a mesh on region [0,1]x[0,1] using quad elements with 9 nodes.
    The numbering for both nodal points and elements is row by row.
    For example a 2x2 mesh is numbered as follows
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

    n_x, n_y = num_el

    numnod = 9 # number of nodes per element
    nn1row = 2 * n_x + 1 # number of nodes in one row
    nn1col = 2 * n_y + 1 # number of nodes in one column

    # create mesh and fill structure
    nelem = n_x * n_y
    nnodes = nn1row * nn1col
    mesh = Mesh(ndim=2, nnodes=nnodes, elshape=6, nelem=nelem, 
                elnumnod=numnod, npoints=4, ncurves=4, 
                topology=np.zeros((numnod, nelem),dtype=int),
                coor=np.zeros((nnodes, 2)), points=np.zeros((4),dtype=int))

    # topology
    for i in range(n_x):
        for j in range(n_y):
            nn1 = 2 * j * nn1row + 2 * i # number of nodes before element (i,j)
            nn2 = nn1 + nn1row # number of nodes before element (i,j) + 1 row
            nn3 = nn2 + nn1row # number of nodes before element (i,j) + 2 rows

            elem = i + j * n_x # element number

            mesh.topology[:, elem] = [nn1, nn1+1, nn1+2, nn2+2, 
                                      nn3+2, nn3+1, nn3, nn2, nn2+1]

    # coordinates
    if ratio is None:
        # equidistant
        deltax = 0.5 / n_x
        deltay = 0.5 / n_y
        for i in range(nn1row):
            for j in range(nn1col):
                node = i + j * nn1row
                mesh.coor[node] = [i * deltax, j * deltay]
    else:
        # non-equidistant
        x1, x3 = np.zeros(nn1row), np.zeros(nn1row),
        x2, x4 = np.zeros(nn1col), np.zeros(nn1col)
        arrays = [x1, x2, x3, x4]
        for idx, array in enumerate(arrays):
            array[::2] = distribute_elements(n_x if idx % 2 == 0 
                                             else n_y, ratio[idx], factor[idx])
        
        # mid-side nodes
        for elem in range(n_x):
            k = 2 * elem # nodes before element elem
            x1[k + 1] = 0.5 * (x1[k] + x1[k + 2])
            x3[k + 1] = 0.5 * (x3[k] + x3[k + 2])
        
        for elem in range(n_y):
            k = 2 * elem # nodes before element elem
            x2[k + 1] = 0.5 * (x2[k] + x2[k + 2])
            x4[k + 1] = 0.5 * (x4[k] + x4[k + 2])

        x3 = 1 - x3[::-1]
        x4 = 1 - x4[::-1]

        # create straight lines in reference square [0,1]x[0,1]
        i_vals = np.arange(nn1row)
        j_vals = np.arange(nn1col)

        I, J = np.meshgrid(i_vals, j_vals)

        nodes = I + J * nn1row
        D = (x1[I] - x3[I]) * (x4[J] - x2[J]) - 1

        mesh.coor[nodes, 0] = (x4[J] * (x1[I] - x3[I]) - x1[I]) / D
        mesh.coor[nodes, 1] = (x1[I] * (x4[J] - x2[J]) - x4[J]) / D

    # points
    mesh.points = np.array([0, nn1row-1, 
                            nn1row*nn1col-1, nn1row*(nn1col-1)],dtype=int)

    # curves
    temp_a = [nn1row, nn1col, nn1row, nn1col]
    temp_b = [n_x, n_y, n_x, n_y]

    mesh.curves = [Geometry(elshape=2,ndim=2,elnumnod=3,nnodes=temp_a[curve],
                            nelem=temp_b[curve]) for curve in range(4)]

    # nodes of curves
    mesh.curves[0].nodes = np.arange(nn1row)
    mesh.curves[1].nodes = nn1row * ( np.arange(nn1col) + 1 ) - 1
    mesh.curves[2].nodes = (nn1col - 1) * nn1row + np.arange(nn1row)[::-1]
    mesh.curves[3].nodes = nn1row * np.arange(nn1col)[::-1]

    # topology of elements on curves
    for curve in mesh.curves:

        nelem = curve.nelem
        elnumnod = curve.elnumnod
        curve.topology = np.zeros((elnumnod, nelem, 2), dtype=int)

        # local numbering
        for elem in range(nelem):
            curve.topology[:, elem, 0] = np.arange(2 * elem, 2 * elem + 3)

        # global numbering
        curve.topology[:, :, 1] = curve.nodes[ curve.topology[:, :, 0] ]

    return mesh