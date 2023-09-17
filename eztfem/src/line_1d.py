import numpy as np
from eztfem.src.mesh_class import Mesh, Geometry
from eztfem.src.distribute_elements import distribute_elements

def line1d(ne, eltype, **kwargs):
    '''
    LINE1D  Simple meshgenerator for 1D lines on the interval [0:1].
      Simple meshgenerator for 1D lines on the interval [0:1] with an optional
      translation and scaling of the region.
      mesh = LINE1D ( n, eltype, 'option1', value1, .... )
      input:
        ne: number of elements 
        eltype: shape number 
                 eltype='line2': 2 node elements
                 eltype='line3': 3 node elements
      optional arguments:
        string, value couples to set options:
        'origin', origin of the domain
        'length', length of the domain
        'ratio', where
            ratio = 0: equidistant mesh
            ratio = 1: the size of the last element is factor times the first
            ratio = 2: the size of an element is factor times the previous one
            ratio = 3: the size of the last element is 1/factor times the first
            ratio = 4: the size of an element is 1/factor times the previous one
           default=0
        'factor' factor
           default=1
      For example:
        mesh = line1d ( 10, 1, 'origin', 1, 'length', 2 )
      to create a 10-element line mesh with 2-node elements on the domain [1,3].
      output:
        mesh: mesh structure. See function quadrilateral2d for the components.
              Note, that the components curves are not present for LINE1D!
    '''

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
      Nodes:    
          1 --  2 --  3 --  4 --  5    
      Elements:    
          x --  x --  x --  x --  x
             1     2     3     4    
    Also generated are two points and the end of the interval
    
         P1 --------------------- P2
   '''
    
    print('line1d_2node')
    
    # create mesh object
    mesh = Mesh(ndim=1,nnodes=n+1,elshape=1,nelem=n,elnumnod=2,npoints=2,
                topology=np.zeros((2, n)), coor=np.zeros((n + 1, 1)),
                points = np.zeros(2))

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
    Nodes:    
        1 --  2 --  3 --  4 --  5    
    Elements:    
        x --  x --  x --  x --  x
              1           2     
    Also generated are two points and the end of the interval    
       P1 --------------------- P2
    '''
    print('line1d_3node')
    
    # create mesh object
    mesh = Mesh(ndim=1,nnodes=2*n+1,elshape=2,nelem=n,elnumnod=3,npoints=2,
                topology=np.zeros((3, n)), coor=np.zeros((2*n + 1, 1)),
                points = np.zeros(2))

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