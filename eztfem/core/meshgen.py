import numpy as np
from .mesh import Mesh, Geometry


def distribute_elements(nelem, ratio, factor):
    """
    Generate non-equidistant n elements on the interval [0, 1].

    Parameters
    ----------
    n : int
        Number of elements.
    ratio : int
        Determines the distribution of elements:
        0: Equidistant mesh.
        1: The size of the last element is `factor` times the first.
        2: The size of an element is `factor` times the previous one.
        3: The size of the last element is `1/factor` times the first.
        4: The size of an element is `1/factor` times the previous one.
    factor : float
        Factor used in the distribution calculations.

    Returns
    -------
    x : numpy.ndarray
        Coordinates of n+1 points.

    Notes
    -----
    The interval [0, 1] is divided into n elements:

    .. math::
        dx_{i+1} = g \\cdot dx_{i-1}

    .. math::
        1 = (1 + g + g^2 + g^3 + \\ldots + g^{n-1}) \\cdot dx_1 =
        \\frac{1 - g^n}{1 - g} \\cdot dx_1

    With :math:`fac = 1 + g + g^2 + \\ldots + g^{n-1}`, we have:

    .. math::
        dx_1 = \\frac{1}{fac} = \\frac{1 - g}{1 - g^n}

    .. math::
        dx_n = g^{n-1} \\cdot dx_1

    """

    if factor < 0:
        raise ValueError("Negative factor")

    if nelem <= 1:
        raise ValueError("Number of elements must be at least 2")

    match ratio:
        case 0:
            g = 1
        case 1:
            g = np.exp(np.log(factor) / (nelem - 1))
        case 2:
            g = factor
        case 3:
            g = np.exp(-np.log(factor) / (nelem - 1))
        case 4:
            g = 1 / factor
        case _:
            raise ValueError(f"Invalid value for ratio: {ratio}")

    # generate mesh
    fac = 1.0
    for i in range(1, nelem):
        fac = fac + g**i

    # size of first element
    dx_1 = 1.0 / fac

    # generate all elements
    x = np.zeros(nelem+1)
    x[0] = 0
    dx = dx_1
    for i in range(1, nelem+1):
        x[i] = x[i-1] + dx
        dx = g * dx

    # test whether x(n+1) = 1
    if abs(x[-1] - 1) > 1e-10:
        raise ValueError(f"End value x(n+1) != 1: {x[-1]:.5e}")
    else:
        x[-1] = 1.0

    return x


def line1d(ne, eltype, **kwargs):
    """
    Simple mesh generator for 1D lines on the interval [0, 1].

    Generates a simple 1D line mesh on the interval [0, 1] with optional
    translation and scaling of the region.

    Parameters
    ----------
    ne : int
        Number of elements.
    eltype : str
        Shape number.
        - 'line2' : 2-node elements.
        - 'line3' : 3-node elements.

    Keyword arguments
    -----------------
    origin : float, optional
        Origin of the domain.
    length : float, optional
        Length of the domain.
    ratio : int, optional, default=0
        Mesh ratio.
        - 0 : Equidistant mesh.
        - 1 : The size of the last element is factor times the first.
        - 2 : The size of an element is factor times the previous one.
        - 3 : The size of the last element is 1/factor times the first.
        - 4 : The size of an element is 1/factor times the previous one.
    factor : float, optional, default=1
        Factor for the mesh ratio.

    Returns
    -------
    mesh : Mesh
        Mesh object. See function `quadrilateral2d` for the components.
        Note that the components curves are not present for `LINE1D`.

    Examples
    --------
    >>> mesh = line1d(10, 'line2', origin=1, length=2)
    Creates a 10-element line mesh with 2-node elements on the domain [1, 3].

    """

    # optional arguments
    ori = kwargs.get('origin', None)
    length = kwargs.get('length', None)
    ratio = kwargs.get('ratio', None)
    factor = kwargs.get('factor', None)

    # mesh
    if eltype == 'line2':
        mesh = line1d_2node(ne, ratio, factor)
    elif eltype == 'line3':
        mesh = line1d_3node(ne, ratio, factor)
    else:
        raise ValueError(f'Invalid eltype = {eltype}')

    mesh.ncurves = 0

    # translate and scale unit region
    if length is not None:
        mesh.coor[:, 0] *= length
    if ori is not None:
        mesh.coor[:, 0] += ori

    return mesh


def line1d_2node(n, ratio, factor):
    '''
    Generate a mesh on region [0,1] using line elements with 2 nodes.
    The numbering is straightforward:
    For example a 4 element mesh is numbered as follows

    ..

      Nodes:
          1 --  2 --  3 --  4 --  5
      Elements:
          x --  x --  x --  x --  x
             1     2     3     4

    Also generated are two points and the end of the interval

    ..

         P1 --------------------- P2

    '''

    print('line1d_2node')

    # create mesh object
    mesh = Mesh(ndim=1, nnodes=n+1, elshape=1, nelem=n, elnumnod=2, npoints=2,
                topology=np.zeros((2, n)), coor=np.zeros((n + 1, 1)),
                points=np.zeros(2))

    # topology
    for elem in range(n):
        mesh.topology[0, elem] = elem
        mesh.topology[1, elem] = elem + 1

    # coordinates
    if not ratio:
        # equidistant
        deltax = 1 / n
        for node in range(mesh.nnodes):
            mesh.coor[node, 0] = node * deltax
    else:
        # non-equidistant
        x = distribute_elements(n, ratio, factor)
        mesh.coor[:, 0] = x

    # points
    mesh.points[0] = 0
    mesh.points[1] = mesh.nnodes - 1

    return mesh


def line1d_3node(n, ratio, factor):
    '''
    Generate a mesh on region [0,1] using line elements with 3 nodes.
    The numbering is straightforward:
    For example a 2 element mesh is numbered as follows

    ..

    Nodes:
        1 --  2 --  3 --  4 --  5
    Elements:
        x --  x --  x --  x --  x
              1           2

    Also generated are two points and the end of the interval

    ..

       P1 --------------------- P2

    '''

    print('line1d_3node')

    # create mesh object
    mesh = Mesh(ndim=1, nnodes=2*n+1, elshape=2, nelem=n, elnumnod=3,
                npoints=2, topology=np.zeros((3, n)),
                coor=np.zeros((2*n + 1, 1)), points=np.zeros(2))

    # topology
    for elem in range(n):
        mesh.topology[0, elem] = 2*elem
        mesh.topology[1, elem] = 2*elem + 1
        mesh.topology[2, elem] = 2*elem + 2

    # coordinates
    if not ratio:
        # equidistant
        deltax = 0.5 / n
        for node in range(mesh.nnodes):
            mesh.coor[node, 0] = node * deltax
    else:
        # non-equidistant
        x = distribute_elements(n, ratio, factor)
        mesh.coor[0, 0] = x[0]
        for elem in range(mesh.nelem):
            k = 2 * elem
            mesh.coor[k + 1, 0] = (x[elem+1] + x[elem]) / 2
            mesh.coor[k + 2, 0] = x[elem + 1]

    # points
    mesh.points[0] = 0
    mesh.points[1] = mesh.nnodes - 1

    return mesh


def quadrilateral2d(num_el, eltype, **kwargs):
    """
    Simple mesh generator for quadrilateral 2D regions.

    Generates a simple 2D mesh on the region [0, 1] x [0, 1] with optional
    translation and/or scaling or deformation of the region.

    Parameters
    ----------
    num_el : array_like
        Number of elements in x and y direction, e.g., [nx, ny].
    eltype : str
        Element type.
        - 'tria3' : 3-node triangle.
        - 'tria4' : 4-node triangle.
        - 'tria6' : 6-node triangle.
        - 'tria7' : 7-node triangle.
        - 'quad4' : 4-node quadrilateral.
        - 'quad9' : 9-node quadrilateral.
        - 'quad5' : 5-node quadrilateral.

    Keyword arguments
    -----------------
    origin : array_like, optional, default=[0, 0]
        Origin of the domain [ox, uy].
    length : array_like, optional, default=[1.0, 1.0]
        Length and width of the domain [lx, ly].
    vertices : array_like, optional, default=[[0, 0], [1, 0], [1, 1], [0, 1]]
        Four vertices of the domain in the format
        [[x1, y1], [x2, y2], [x3, y3], [x4, y4]].
        Note that origin/length arguments are ignored if vertices are
        given.
    ratio : list of int, optional, default=[0, 0, 0, 0]
        A vector [ratio1, ratio2, ratio3, ratio4] where for each of the
        four curves the ratio is given:
        - 0 : Equidistant mesh.
        - 1 : The size of the last element is factor times the first.
        - 2 : The size of an element is factor times the previous one.
        - 3 : The size of the last element is 1/factor times the first.
        - 4 : The size of an element is 1/factor times the previous one.
    factor : array_like, optional, default=[1, 1, 1, 1]
        A vector [factor1, factor2, factor3, factor4] where for each of
        the four curves the factor is given.

    Returns
    -------
    mesh : Mesh
        Mesh object

    Examples
    --------
    >>> mesh = quadrilateral2d([10, 10], 'tria3', origin=[1, 2], length=[2, 3])
    Creates a 10x10 mesh of 3-node triangles, where the left lower corner of
    the domain is (1, 2) and the size of the domain is 2 by 3.

    >>> mesh = quadrilateral2d([10, 10], 'quad5',
                               vertices=[[1, 1], [3, 2], [4, 5], [-1, 4]])
    Creates a 10x10 mesh of 4-node quadrilateral elements, where the left
    lower corner, right lower corner, upper right corner, and upper left
    corner of the domain are (1, 1), (3, 2), (4, 5), and (-1, 4), respectively.
    The edges are straight.

    """

    # optional arguments
    ori = kwargs.get('origin', None)
    length = kwargs.get('length', None)
    def_verts = kwargs.get('vertices', None)
    ratio = kwargs.get('ratio', None)
    factor = kwargs.get('factor', None)

    # mesh
    if eltype == 'tria3':  # 3 node triangle
        mesh = rectangle2d_tria3(num_el, ratio, factor)
    elif eltype == 'tria4':  # 4 node triangle
        mesh = rectangle2d_tria4(num_el, ratio, factor)
    elif eltype == 'tria6':  # 6 node triangle
        mesh = rectangle2d_tria6(num_el, ratio, factor)
    elif eltype == 'tria7':  # 7 node triangle
        mesh = rectangle2d_tria7(num_el, ratio, factor)
    elif eltype == 'quad4':  # 4 node quadrilateral
        mesh = rectangle2d_quad4(num_el, ratio, factor)
    elif eltype == 'quad5':  # 5 node quadrilateral
        mesh = rectangle2d_quad5(num_el, ratio, factor)
    elif eltype == 'quad9':  # 9 node quadrilateral
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
        coor_x_1 = x_weights_1 * def_verts[0, :] \
            + x_weights_2 * def_verts[1, :]
        coor_x_2 = x_weights_1 * def_verts[3, :] \
            + x_weights_2 * def_verts[2, :]
        mesh.coor = y_weights_1 * coor_x_1 + y_weights_2 * coor_x_2
    else:
        if length is not None:
            mesh.coor *= length
        if ori is not None:
            mesh.coor += ori

    return mesh


def rectangle2d_tria3(num_el, ratio, factor):
    """
    Generate a mesh on region [0,1]x[0,1] using triangular elements with 3
    nodes.

    The numbering for both nodal points and elements is row by row.
    For example a 2x2 mesh is numbered as follows

    ::

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

    ::

                      <-----
                         C3
             P4 --------------------- P3
              |                       |           Note:   C1=P1-P2
              |                       |                   C2=P2-P3
         |    |                       |    ^              C3=P3-P4
         | C4 |                       | C2 |              C4=P4-P1
        \\|/   |                       |    |
              |                       |
              |                       |
             P1 --------------------- P2
                         C1
                       ----->

    """

    print('rectangle2d_tria3')

    n_x, n_y = num_el

    numnod = 3  # number of nodes per element
    nn1row = n_x + 1  # number of nodes in one row
    nn1col = n_y + 1  # number of nodes in one column

    # create mesh and fill structure
    nelem = 2 * n_x * n_y
    nnodes = nn1row * nn1col
    mesh = Mesh(ndim=2, nnodes=nnodes, elshape=3, nelem=nelem,
                elnumnod=numnod, npoints=4, ncurves=4,
                topology=np.zeros((numnod, nelem), dtype=int),
                coor=np.zeros((nnodes, 2)), points=np.zeros((4), dtype=int))

    # topology
    for i in range(n_x):
        for j in range(n_y):
            nn1 = j * nn1row + i  # number of nodes before element (i,j)
            nn2 = nn1 + nn1row  # number of nodes before element (i,j) + 1 row

            elem = 2 * i + 2 * j * n_x  # element number
            mesh.topology[:, elem] = [nn1, nn2+1, nn2]
            elem += 1  # element number
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
            array[:] = distribute_elements(n_x
                                           if idx % 2 == 0
                                           else n_y, ratio[idx], factor[idx])

        x3 = 1 - x3[::-1]
        x4 = 1 - x4[::-1]

        # create straight lines in reference square [0,1]x[0,1]
        for i in range(nn1row):
            for j in range(nn1col):
                node = i + j*nn1row
                D = (x1[i] - x3[i]) * (x4[j] - x2[j]) - 1
                mesh.coor[node, 0] = (x4[j] * (x1[i] - x3[i]) - x1[i]) / D
                mesh.coor[node, 1] = (x1[i] * (x4[j] - x2[j]) - x4[j]) / D

    # points
    mesh.points = np.array([0, nn1row-1,
                            nn1row*nn1col-1, nn1row*(nn1col-1)], dtype=int)

    # curves
    temp_a = [nn1row, nn1col, nn1row, nn1col]
    temp_b = [n_x, n_y, n_x, n_y]

    mesh.curves = [Geometry(elshape=1, ndim=2, elnumnod=2,
                            nnodes=temp_a[curve], nelem=temp_b[curve])
                   for curve in range(4)]

    # nodes of curves
    mesh.curves[0].nodes = np.arange(nn1row)
    mesh.curves[1].nodes = nn1row * (np.arange(nn1col) + 1) - 1
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
        curve.topology[:, :, 1] = curve.nodes[curve.topology[:, :, 0]]

    return mesh


def rectangle2d_tria4(num_el, ratio, factor):
    """
    Generate a mesh on region [0,1]x[0,1] using triangular elements with 4
    nodes.

    The numbering for both nodal points and elements is row by row.
    For example a 2x2 mesh is numbered as follows

    ::

      Nodes:

          15------- 16------- 17
          | 11    / | 13    / |
          |    /    |    /    |
          | /    12 | /    14 |
          8 ------- 9 ------- 10
          |  4    / | 6     / |
          |    /    |    /    |
          | /    5  | /     7 |
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

    ::

                      <-----
                         C3
             P4 --------------------- P3
              |                       |           Note:   C1=P1-P2
              |                       |                   C2=P2-P3
         |    |                       |    ^              C3=P3-P4
         | C4 |                       | C2 |              C4=P4-P1
        \\|/   |                       |    |
              |                       |
              |                       |
             P1 --------------------- P2
                         C1
                       ----->

    """

    print('rectangle2d_tria4')

    n_x, n_y = num_el

    numnod = 4  # number of nodes per element
    nn1row = n_x + 1  # number of nodes in one row
    nn1col = n_y + 1  # number of nodes in one column

    # create mesh and fill structure
    nelem = 2 * n_x * n_y
    nnodes = nn1row * nn1col + 2 * n_x * n_y
    mesh = Mesh(ndim=2, nnodes=nnodes, elshape=10, nelem=nelem,
                elnumnod=numnod, npoints=4, ncurves=4,
                topology=np.zeros((numnod, nelem), dtype=int),
                coor=np.zeros((nnodes, 2)), points=np.zeros((4), dtype=int))

    # topology
    for i in range(n_x):
        for j in range(n_y):
            # number of nodes before element (i,j)
            nn1 = j * (nn1row+2*n_x) + i

            # number of nodes before element (i,j) + 1 row
            nn2 = nn1 + nn1row + i

            # number of nodes before element (i,j) + 2 rows
            nn3 = nn1 + nn1row + 2*n_x

            elem = 2 * i + 2 * j * n_x  # element number
            mesh.topology[:, elem] = [nn1, nn3+1, nn3, nn2]
            elem += 1  # element number
            mesh.topology[:, elem] = [nn1, nn1+1, nn3+1, nn2+1]

    # coordinates
    if ratio is None:
        # equidistant
        deltax = 1.0 / n_x
        deltay = 1.0 / n_y
        for i in range(nn1row):
            for j in range(nn1col):
                node = i + j * (nn1row + 2 * n_x)
                mesh.coor[node] = [i * deltax, j * deltay]
        for i in range(n_x):
            for j in range(n_y):
                node = 2*i + j*(nn1row+2*n_x) + nn1row
                mesh.coor[node, 0] = i * deltax + deltax/3
                mesh.coor[node, 1] = j * deltay + 2*deltay/3
                node = node + 1
                mesh.coor[node, 0] = i * deltax + 2*deltax/3
                mesh.coor[node, 1] = j * deltay + deltay/3
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
        for i in range(nn1row):
            for j in range(nn1col):
                node = i + j*(nn1row+2*n_x)
                D = (x1[i] - x3[i]) * (x4[j] - x2[j]) - 1
                mesh.coor[node, 0] = (x4[j] * (x1[i] - x3[i]) - x1[i]) / D
                mesh.coor[node, 1] = (x1[i] * (x4[j] - x2[j]) - x4[j]) / D

        for i in range(n_x):
            for j in range(n_y):
                node = 2*i + j*(nn1row+2*n_x) + nn1row
                node1 = i + j*(nn1row+2*n_x)
                node3 = i + (j+1)*(nn1row+2*n_x)
                mesh.coor[node, :] = (mesh.coor[node1, :]
                                      + mesh.coor[node3, :]
                                      + mesh.coor[node3+1, :]) / 3
                node = node + 1
                node3 = node3 + 1
                mesh.coor[node, :] = (mesh.coor[node1, :]
                                      + mesh.coor[node1+1, :]
                                      + mesh.coor[node3, :]) / 3

    # points
    mesh.points = np.array([0, nn1row-1, nn1row*nn1col+2*n_x*n_y-1,
                            nn1row*(nn1col-1)+2*n_x*n_y], dtype=int)

    # curves
    temp_a = [nn1row, nn1col, nn1row, nn1col]
    temp_b = [n_x, n_y, n_x, n_y]

    mesh.curves = [Geometry(elshape=1, ndim=2, elnumnod=2,
                            nnodes=temp_a[curve], nelem=temp_b[curve])
                   for curve in range(4)]

    # nodes of curves
    mesh.curves[0].nodes = np.arange(nn1row)
    mesh.curves[1].nodes = (nn1row+2*n_x) * (np.arange(nn1col) + 1) - 2*n_x - 1
    mesh.curves[2].nodes = (nn1col - 1) * nn1row + 2*n_x*n_y \
        + np.arange(nn1row)[::-1]
    mesh.curves[3].nodes = (nn1row+2*n_x) * np.arange(nn1col)[::-1]

    # topology of elements on curves
    for curve in mesh.curves:

        nelem = curve.nelem
        elnumnod = curve.elnumnod
        curve.topology = np.zeros((elnumnod, nelem, 2), dtype=int)

        # local numbering
        for elem in range(nelem):
            curve.topology[:, elem, 0] = np.arange(elem, elem + 2)

        # global numbering
        curve.topology[:, :, 1] = curve.nodes[curve.topology[:, :, 0]]

    return mesh


def rectangle2d_tria6(num_el, ratio, factor):
    """
    Generate a mesh on region [0,1]x[0,1] using triangular elements with 6
    nodes.

    The numbering for both nodal points and elements is row by row.
    For example a 2x2 mesh is numbered as follows

    ::

      Nodes:

         21 -- 22 -- 23 -- 24 -- 25
          |        /  |        /  |
         16    17    18    19    20
          | /         | /         |
         11 -- 12 -- 13 -- 14 -- 15
          |        /  |        /  |
          6     7     8     9    10
          | /         | /         |
          1 --  2 --  3 --  4 --  5

      Elements:

          x --  x --  x --  x --  x
          | 5      /  | 7      /  |
          x     x     x     x     x
          | /      6  | /      8  |
          x --  x --  x --  x --  x
          |  1     /  |  3     /  |
          x     x     x     x     x
          | /      2  | /      4  |
          x --  x --  x --  x --  x

    Also generated are four points and four curves:

    ::

                      <-----
                         C3
             P4 --------------------- P3
              |                       |           Note:   C1=P1-P2
              |                       |                   C2=P2-P3
         |    |                       |    ^              C3=P3-P4
         | C4 |                       | C2 |              C4=P4-P1
        \\|/   |                       |    |
              |                       |
              |                       |
             P1 --------------------- P2
                         C1
                       ----->

    """

    print('rectangle_tria6')

    n_x, n_y = num_el

    numnod = 6  # number of nodes per element
    nn1row = 2 * n_x + 1  # number of nodes in one row
    nn1col = 2 * n_y + 1  # number of nodes in one column

    # create mesh and fill structure
    nelem = 2 * n_x * n_y
    nnodes = nn1row * nn1col
    mesh = Mesh(ndim=2, nnodes=nnodes, elshape=4, nelem=nelem,
                elnumnod=numnod, npoints=4, ncurves=4,
                topology=np.zeros((numnod, nelem), dtype=int),
                coor=np.zeros((nnodes, 2)), points=np.zeros((4), dtype=int))

    # topology
    for i in range(n_x):
        for j in range(n_y):
            # number of nodes before element (i,j)
            nn1 = 2 * j * nn1row + 2 * i

            nn2 = nn1 + nn1row  # number of nodes before element (i,j) + 1 row
            nn3 = nn2 + nn1row  # number of nodes before element (i,j) + 2 rows

            elem = 2 * i + 2 * j * n_x  # element number

            mesh.topology[:, elem] = [nn1, nn2+1, nn3+2, nn3+1, nn3, nn2]

            elem += 1  # element number

            mesh.topology[:, elem] = [nn1, nn1+1, nn1+2, nn2+2, nn3+2, nn2+1]

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
            k = 2 * elem  # nodes before element elem
            x1[k + 1] = 0.5 * (x1[k] + x1[k + 2])
            x3[k + 1] = 0.5 * (x3[k] + x3[k + 2])

        for elem in range(n_y):
            k = 2 * elem  # nodes before element elem
            x2[k + 1] = 0.5 * (x2[k] + x2[k + 2])
            x4[k + 1] = 0.5 * (x4[k] + x4[k + 2])

        x3 = 1 - x3[::-1]
        x4 = 1 - x4[::-1]

        # create straight lines in reference square [0,1]x[0,1]
        for i in range(nn1row):
            for j in range(nn1col):
                node = i + j*nn1row
                D = (x1[i] - x3[i]) * (x4[j] - x2[j]) - 1
                mesh.coor[node, 0] = (x4[j] * (x1[i] - x3[i]) - x1[i]) / D
                mesh.coor[node, 1] = (x1[i] * (x4[j] - x2[j]) - x4[j]) / D

    # points
    mesh.points = np.array([0, nn1row-1,
                            nn1row*nn1col-1, nn1row*(nn1col-1)], dtype=int)

    # curves
    temp_a = [nn1row, nn1col, nn1row, nn1col]
    temp_b = [n_x, n_y, n_x, n_y]

    mesh.curves = [Geometry(elshape=2, ndim=2, elnumnod=3,
                            nnodes=temp_a[curve],
                            nelem=temp_b[curve]) for curve in range(4)]

    # nodes of curves
    mesh.curves[0].nodes = np.arange(nn1row)
    mesh.curves[1].nodes = nn1row * (np.arange(nn1col) + 1) - 1
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
        curve.topology[:, :, 1] = curve.nodes[curve.topology[:, :, 0]]

    return mesh


def rectangle2d_tria7(num_el, ratio, factor):
    """
    Generate a mesh on region [0,1]x[0,1] using triangular elements with 7
    nodes.

    The numbering for both nodal points and elements is row by row.
    For example a 2x2 mesh is numbered as follows

    ::

      Nodes:

         29 -- 30 -- 31 -- 32 -- 25
          | 20     /  | 22     /  |
         24    25    26    27    28
          | /     21  | /     23  |
         15 -- 16 -- 17 -- 18 -- 19
          | 6      /  | 8      /  |
          10   11    12    13    14
          | /      7  | /      9  |
          1 --  2 --  3 --  4 --  5

      Elements:

          x --  x --  x --  x --  x
          | 5      /  | 7      /  |
          x     x     x     x     x
          | /      6  | /      8  |
          x --  x --  x --  x --  x
          |  1     /  |  3     /  |
          x     x     x     x     x
          | /      2  | /      4  |
          x --  x --  x --  x --  x

    Also generated are four points and four curves:

    ::

                      <-----
                         C3
             P4 --------------------- P3
              |                       |           Note:   C1=P1-P2
              |                       |                   C2=P2-P3
         |    |                       |    ^              C3=P3-P4
         | C4 |                       | C2 |              C4=P4-P1
        \\|/   |                       |    |
              |                       |
              |                       |
             P1 --------------------- P2
                         C1
                       ----->

    """

    print('rectangle_tria7')

    n_x, n_y = num_el

    numnod = 7  # number of nodes per element
    nn1row = 2 * n_x + 1  # number of nodes in one row
    nn1col = 2 * n_y + 1  # number of nodes in one column

    # create mesh and fill structure
    nelem = 2 * n_x * n_y
    nnodes = nn1row * nn1col + 2*n_x*n_y
    mesh = Mesh(ndim=2, nnodes=nnodes, elshape=7, nelem=nelem,
                elnumnod=numnod, npoints=4, ncurves=4,
                topology=np.zeros((numnod, nelem), dtype=int),
                coor=np.zeros((nnodes, 2)), points=np.zeros((4), dtype=int))

    # topology
    for i in range(n_x):
        for j in range(n_y):
            # number of nodes before element (i,j)
            nn1 = j * (2*nn1row + 2*n_x) + 2*i
            # number of nodes before element (i,j) + 1 row
            nn2 = nn1 + nn1row
            # number of nodes before element (i,j) + 2 rows
            nn3 = nn2 + 2*n_x
            # number of nodes before element (i,j) + 3 rows
            nn4 = nn3 + nn1row

            elem = 2*i + j*2*n_x   # element number

            mesh.topology[:, elem] = [nn1, nn3+1, nn4+2, nn4+1, nn4, nn3, nn2]

            elem += 1  # element number

            mesh.topology[:, elem] = [nn1, nn1+1, nn1+2, nn3+2, nn4+2, nn3+1,
                                      nn2+1]

    # coordinates
    if ratio is None:
        # equidistant
        deltax = 1.0 / n_x
        deltay = 1.0 / n_y
        for i in range(nn1row):
            for j in range(n_y+1):
                node = i + j*(2*nn1row+2*n_x)
                mesh.coor[node, 0] = i * deltax/2
                mesh.coor[node, 1] = j * deltay
        for i in range(n_x):
            for j in range(n_y):
                node = 2*i + nn1row + j*(2*nn1row+2*n_x)
                mesh.coor[node, 0] = i * deltax + deltax/3
                mesh.coor[node, 1] = j * deltay + 2*deltay/3
                node = node + 1
                mesh.coor[node, 0] = i * deltax + 2*deltax/3
                mesh.coor[node, 1] = j * deltay + deltay/3
        for i in range(nn1row):
            for j in range(n_y):
                node = i + nn1row+2*n_x + j*(2*nn1row+2*n_x)
                mesh.coor[node, 0] = i * deltax/2
                mesh.coor[node, 1] = j * deltay + deltay/2
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
            k = 2 * elem  # nodes before element elem
            x1[k + 1] = 0.5 * (x1[k] + x1[k + 2])
            x3[k + 1] = 0.5 * (x3[k] + x3[k + 2])

        for elem in range(n_y):
            k = 2 * elem  # nodes before element elem
            x2[k + 1] = 0.5 * (x2[k] + x2[k + 2])
            x4[k + 1] = 0.5 * (x4[k] + x4[k + 2])

        x3 = 1 - x3[::-1]
        x4 = 1 - x4[::-1]

        # create straight lines in reference square [0,1]x[0,1]
        for i in range(nn1row):
            for j in range(nn1col):
                if j % 2 == 1:  # even row (Python 0-based indexing)
                    node = i + j*nn1row + (j+1)*n_x
                else:  # odd row
                    node = i + j*(nn1row+n_x)
                D = (x1[i] - x3[i]) * (x4[j] - x2[j]) - 1
                mesh.coor[node, 0] = (x4[j] * (x1[i] - x3[i]) - x1[i]) / D
                mesh.coor[node, 1] = (x1[i] * (x4[j] - x2[j]) - x4[j]) / D

        for i in range(n_x):
            for j in range(n_y):
                node = 2*i + nn1row + j*(2*nn1row+2*n_x)
                node1 = 2*i + j*(2*nn1row+2*n_x)
                node3 = 2*i + (j+1)*(2*nn1row+2*n_x)
                mesh.coor[node, :] = (mesh.coor[node1, :]
                                      + mesh.coor[node3, :]
                                      + mesh.coor[node3+2, :]) / 3
                node = node + 1
                node3 = node3 + 2
                mesh.coor[node, :] = (mesh.coor[node1, :]
                                      + mesh.coor[node1+2, :]
                                      + mesh.coor[node3, :]) / 3

    # points
    mesh.points = np.array([0, nn1row-1,  nn1row*nn1col+2*n_x*n_y-1,
                            nn1row*(nn1col-1)+2*n_x*n_y], dtype=int)

    # curves
    temp_a = [nn1row, nn1col, nn1row, nn1col]
    temp_b = [n_x, n_y, n_x, n_y]

    mesh.curves = [Geometry(elshape=2, ndim=2, elnumnod=3,
                            nnodes=temp_a[curve],
                            nelem=temp_b[curve]) for curve in range(4)]

    # nodes of curves

    # curve 1
    mesh.curves[0].nodes = np.arange(nn1row)

    # curve 2
    mesh.curves[1].nodes = np.zeros(nn1col, dtype=int)
    mesh.curves[1].nodes[0:nn1col:2] = \
        (nn1row + n_x) * np.arange(1, nn1col+1, 2) - n_x - 1
    mesh.curves[1].nodes[1:nn1col:2] = \
        (nn1row + n_x) * np.arange(2, nn1col+1, 2) - 1

    # curve 3
    mesh.curves[2].nodes = \
        (nn1col-1)*nn1row + 2*n_x*n_y + np.arange(nn1row, 0, -1) - 1

    # curve 4
    mesh.curves[3].nodes = np.zeros(nn1col, dtype=int)
    mesh.curves[3].nodes[0:nn1col:2] = \
        1-nn1row-n_x + (nn1row+n_x)*(np.arange(nn1col, 0, -2)) - 1
    mesh.curves[3].nodes[1:nn1col:2] = \
        1-nn1row + (nn1row+n_x)*(np.arange(nn1col-1, 0, -2)) - 1

    # topology of elements on curves
    for curve in mesh.curves:

        nelem = curve.nelem
        elnumnod = curve.elnumnod
        curve.topology = np.zeros((elnumnod, nelem, 2), dtype=int)

        # local numbering
        for elem in range(nelem):
            curve.topology[:, elem, 0] = np.arange(2 * elem, 2 * elem + 3)

        # global numbering
        curve.topology[:, :, 1] = curve.nodes[curve.topology[:, :, 0]]

    return mesh


def rectangle2d_quad4(num_el, ratio, factor):
    """
    Generate a mesh on region [0,1]x[0,1] using quad elements with 4 nodes.

    The numbering for both nodal points and elements is row by row.
    For example a 2x2 mesh is numbered as follows

    ::

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

    ::

                      <-----
                         C3
             P4 --------------------- P3
              |                       |           Note:   C1=P1-P2
              |                       |                   C2=P2-P3
         |    |                       |    ^              C3=P3-P4
         | C4 |                       | C2 |              C4=P4-P1
        \\|/   |                       |    |
              |                       |
              |                       |
             P1 --------------------- P2
                         C1
                       ----->

    """

    print('rectangle_quad4')

    n_x, n_y = num_el

    numnod = 4  # number of nodes per element
    nn1row = n_x + 1  # number of nodes in one row
    nn1col = n_y + 1  # number of nodes in one column

    # create mesh and fill structure
    nelem = n_x * n_y
    nnodes = nn1row * nn1col
    mesh = Mesh(ndim=2, nnodes=nnodes, elshape=5, nelem=nelem,
                elnumnod=numnod, npoints=4, ncurves=4,
                topology=np.zeros((numnod, nelem), dtype=int),
                coor=np.zeros((nnodes, 2)), points=np.zeros((4), dtype=int))

    # topology
    for i in range(n_x):
        for j in range(n_y):
            nn1 = j * nn1row + i  # number of nodes before element (i,j)
            nn2 = nn1 + nn1row  # number of nodes before element (i,j) + 1 row

            elem = i + j * n_x  # element number

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
        for i in range(nn1row):
            for j in range(nn1col):
                node = i + j*nn1row
                D = (x1[i] - x3[i]) * (x4[j] - x2[j]) - 1
                mesh.coor[node, 0] = (x4[j] * (x1[i] - x3[i]) - x1[i]) / D
                mesh.coor[node, 1] = (x1[i] * (x4[j] - x2[j]) - x4[j]) / D

    # points
    mesh.points = np.array([0, nn1row-1,
                            nn1row*nn1col-1, nn1row*(nn1col-1)], dtype=int)

    # curves
    temp_a = [nn1row, nn1col, nn1row, nn1col]
    temp_b = [n_x, n_y, n_x, n_y]

    mesh.curves = [Geometry(elshape=1, ndim=2, elnumnod=2,
                            nnodes=temp_a[curve], nelem=temp_b[curve])
                   for curve in range(4)]

    # nodes of curves
    mesh.curves[0].nodes = np.arange(nn1row)
    mesh.curves[1].nodes = nn1row * (np.arange(nn1col) + 1) - 1
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
        curve.topology[:, :, 1] = curve.nodes[curve.topology[:, :, 0]]

    return mesh


def rectangle2d_quad5(num_el, ratio, factor):
    """
    Generate a mesh on region [0,1]x[0,1] using quad elements with 5 nodes.

    The numbering for both nodal points and elements is row by row.
    For example a 2x2 mesh is numbered as follows

    ::

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

    ::

                      <-----
                         C3
             P4 --------------------- P3
              |                       |           Note:   C1=P1-P2
              |                       |                   C2=P2-P3
         |    |                       |    ^              C3=P3-P4
         | C4 |                       | C2 |              C4=P4-P1
        \\|/   |                       |    |
              |                       |
              |                       |
             P1 --------------------- P2
                         C1
                       ----->

    """

    print('rectangle_quad5')

    n_x, n_y = num_el

    numnod = 5  # number of nodes per element
    nn1row = n_x + 1  # number of nodes in one row
    nn1col = n_y + 1  # number of nodes in one column

    # create mesh and fill structure
    nelem = n_x * n_y
    nnodes = nn1row * nn1col + n_x * n_y
    mesh = Mesh(ndim=2, nnodes=nnodes, elshape=9, nelem=nelem,
                elnumnod=numnod, npoints=4, ncurves=4,
                topology=np.zeros((numnod, nelem), dtype=int),
                coor=np.zeros((nnodes, 2)), points=np.zeros((4), dtype=int))

    # topology
    for i in range(n_x):
        for j in range(n_y):
            # number of nodes before element (i,j)
            nn1 = j * (nn1row + n_x) + i

            nn2 = nn1 + nn1row  # number of nodes before element (i,j) + 1 row
            nn3 = nn2 + n_x  # number of nodes before element (i,j) + 2 rows

            elem = i + j * n_x  # element number

            mesh.topology[:, elem] = [nn1, nn1+1, nn3+1, nn3, nn2]

    # coordinates
    if ratio is None:
        # equidistant
        deltax = 1.0 / n_x
        deltay = 1.0 / n_y
        for i in range(nn1row):
            for j in range(nn1col):
                node = i + j * (nn1row + n_x)
                mesh.coor[node] = [i * deltax, j * deltay]
        for i in range(n_x):
            for j in range(n_y):
                node = i + j * (nn1row + n_x) + nn1row
                mesh.coor[node] = [i * deltax + deltax/2,
                                   j * deltay + deltay/2]
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
        for i in range(nn1row):
            for j in range(nn1col):
                node = i + j*(nn1row+n_x)
                D = (x1[i] - x3[i]) * (x4[j] - x2[j]) - 1
                mesh.coor[node, 0] = (x4[j] * (x1[i] - x3[i]) - x1[i]) / D
                mesh.coor[node, 1] = (x1[i] * (x4[j] - x2[j]) - x4[j]) / D
        for i in range(n_x):
            for j in range(n_y):
                node = i + j*(nn1row+n_x) + nn1row
                node1 = i + j*(nn1row+n_x)
                node4 = i + (j+1)*(nn1row+n_x)
                mesh.coor[node, :] = (mesh.coor[node1, :]
                                      + mesh.coor[node1+1, :]
                                      + mesh.coor[node4, :]
                                      + mesh.coor[node4+1, :]) / 4

    # points
    mesh.points = np.array([0, nn1row-1, nn1row*nn1col+n_x*n_y-1,
                            nn1row*(nn1col-1)+n_x*n_y], dtype=int)

    # curves
    temp_a = [nn1row, nn1col, nn1row, nn1col]
    temp_b = [n_x, n_y, n_x, n_y]

    mesh.curves = [Geometry(elshape=1, ndim=2, elnumnod=2,
                            nnodes=temp_a[curve], nelem=temp_b[curve])
                   for curve in range(4)]

    # nodes of curves
    mesh.curves[0].nodes = np.arange(1, nn1row+1) - 1
    mesh.curves[1].nodes = (nn1row + n_x) * (np.arange(1, nn1col+1)) - n_x - 1
    mesh.curves[2].nodes = (nn1col - 1) * nn1row + n_x * n_y \
        + np.arange(1, nn1row+1)[::-1] - 1
    mesh.curves[3].nodes = 1 - nn1row - n_x + (nn1row + n_x) \
        * np.arange(1, nn1col+1)[::-1] - 1

    # topology of elements on curves
    for curve in mesh.curves:

        nelem = curve.nelem
        elnumnod = curve.elnumnod
        curve.topology = np.zeros((elnumnod, nelem, 2), dtype=int)

        # local numbering
        for elem in range(nelem):
            curve.topology[:, elem, 0] = np.arange(elem, elem + 2)

        # global numbering
        curve.topology[:, :, 1] = curve.nodes[curve.topology[:, :, 0]]

    return mesh


def rectangle2d_quad9(num_el, ratio, factor):
    """
    Generate a mesh on region [0,1]x[0,1] using quad elements with 9 nodes.

    The numbering for both nodal points and elements is row by row.
    For example a 2x2 mesh is numbered as follows

    ::

      Nodes
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

    ::

                    <-----
                        C3
            P4 --------------------- P3
            |                       |           Note:   C1=P1-P2
            |                       |                   C2=P2-P3
        |   |                       |    ^              C3=P3-P4
        | C4|                       | C2 |              C4=P4-P1
       \\|/  |                       |    |
            |                       |
            |                       |
            P1 --------------------- P2
                        C1
                    ----->

    """

    print('rectangle_quad9')

    n_x, n_y = num_el

    numnod = 9  # number of nodes per element
    nn1row = 2 * n_x + 1  # number of nodes in one row
    nn1col = 2 * n_y + 1  # number of nodes in one column

    # create mesh and fill structure
    nelem = n_x * n_y
    nnodes = nn1row * nn1col
    mesh = Mesh(ndim=2, nnodes=nnodes, elshape=6, nelem=nelem,
                elnumnod=numnod, npoints=4, ncurves=4,
                topology=np.zeros((numnod, nelem), dtype=int),
                coor=np.zeros((nnodes, 2)), points=np.zeros((4), dtype=int))

    # topology
    for i in range(n_x):
        for j in range(n_y):
            # number of nodes before element (i,j)
            nn1 = 2 * j * nn1row + 2 * i

            nn2 = nn1 + nn1row  # number of nodes before element (i,j) + 1 row
            nn3 = nn2 + nn1row  # number of nodes before element (i,j) + 2 rows

            elem = i + j * n_x  # element number

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
            k = 2 * elem  # nodes before element elem
            x1[k + 1] = 0.5 * (x1[k] + x1[k + 2])
            x3[k + 1] = 0.5 * (x3[k] + x3[k + 2])

        for elem in range(n_y):
            k = 2 * elem  # nodes before element elem
            x2[k + 1] = 0.5 * (x2[k] + x2[k + 2])
            x4[k + 1] = 0.5 * (x4[k] + x4[k + 2])

        x3 = 1 - x3[::-1]
        x4 = 1 - x4[::-1]

        # create straight lines in reference square [0,1]x[0,1]
        for i in range(nn1row):
            for j in range(nn1col):
                node = i + j*nn1row
                D = (x1[i] - x3[i]) * (x4[j] - x2[j]) - 1
                mesh.coor[node, 0] = (x4[j] * (x1[i] - x3[i]) - x1[i]) / D
                mesh.coor[node, 1] = (x1[i] * (x4[j] - x2[j]) - x4[j]) / D

    # points
    mesh.points = np.array([0, nn1row-1,
                            nn1row*nn1col-1, nn1row*(nn1col-1)], dtype=int)

    # curves
    temp_a = [nn1row, nn1col, nn1row, nn1col]
    temp_b = [n_x, n_y, n_x, n_y]

    mesh.curves = [Geometry(elshape=2, ndim=2, elnumnod=3,
                            nnodes=temp_a[curve], nelem=temp_b[curve])
                   for curve in range(4)]

    # nodes of curves
    mesh.curves[0].nodes = np.arange(nn1row)
    mesh.curves[1].nodes = nn1row * (np.arange(nn1col) + 1) - 1
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
        curve.topology[:, :, 1] = curve.nodes[curve.topology[:, :, 0]]

    return mesh
