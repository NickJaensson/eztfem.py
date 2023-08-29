import numpy as np

# class definition for Mesh objects
class Mesh:

    # initialize attributes (will be filled elsewhere)
    def __init__(self, ndim=0, nnodes=0, nelem=0, elshape=0, elnumnod=0, npoints=0, 
                 ncurves=0, topology=None, coor=None, points=None, curves=None):
        
        # it is not recommended to use mutable objects as default values in function def
        # see: https://docs.python.org/3/tutorial/controlflow.html#default-argument-values
        if topology is None:
            topology = np.array([])
        if coor is None:
            coor = np.array([])
        if points is None:
            points = np.array([])
        if curves is None:
            curves = []
        
        self.ndim = ndim
        self.nnodes = nnodes
        self.nelem = nelem
        self.elshape = elshape
        self.elnumnod = elnumnod
        self.npoints = npoints
        self.ncurves = ncurves

        self.topology = topology
        self.coor = coor
        self.points = points
        self.curves = curves

    # equivalence check for testing against Matlab code
    # NOTE: see the matlab file create_pytfem_tests.m for the use of np.squeeze
    def __eq__(self, other):
       check = [self.ndim == other.ndim,
                self.nnodes == other.nnodes,
                np.allclose(np.squeeze(self.coor),other.coor,atol=1e-15,rtol=0),
                self.nelem == other.nelem,
                self.elshape == other.elshape,
                self.elnumnod == other.elnumnod,
                (np.squeeze(self.topology) == other.topology).all(), 
                self.npoints == other.npoints,
                (self.points == other.points).all(),
                self.ncurves == other.ncurves]
       check2 = [self.curves[i]==other.curves[i] for i in range(self.ncurves)]

       # print a warning when not equivalent (for debugging purposes)
       if not ( all(check) and all(check2) ):
           print("WARNING: Meshes not equivalent:")
           print(check)
           print(check2)

       return all(check) and all(check2)

# class definition for Geometry objects
class Geometry:

    # initialize attributes (will be filled elsewhere)
    def __init__(self, elshape=0, ndim=0, elnumnod=0, nnodes=0, nelem=0, 
                 topology=None, nodes=None):
        
        if topology is None:
            topology = np.array([])
        if nodes is None:
            nodes = np.array([])

        self.elshape = elshape
        self.ndim = ndim
        self.elnumnod = elnumnod
        self.nnodes = nnodes
        self.nelem = nelem
        self.topology = topology
        self.nodes = nodes

    # equivalence check for testing against Matlab code
    # NOTE: see the matlab file create_pytfem_tests.m for the use of np.squeeze
    def __eq__(self, other):
        check = [self.ndim == other.ndim,
                 self.elshape == other.elshape,
                 self.elnumnod == other.elnumnod,
                 self.nnodes == other.nnodes,
                 self.nelem == other.nelem,
                 (self.topology == other.topology).all(),
                 (self.nodes == other.nodes).all()]
        
       # print a warning when not equivalent (for debugging purposes)
        if not all(check):
           print("WARNING: Geometries not equivalent:")
           print(check)

        return all(check)
