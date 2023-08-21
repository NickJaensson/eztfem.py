import numpy as np

# class definition for Mesh objects
class Mesh:

    # initialize attributes (will be filled elsewhere)
    def __init__(self):
        self.ndim = 0
        self.nnodes = 0
        self.coor = np.array([])
        self.nelem = 0
        self.elshape = 0
        self.elnumnod = 0
        self.topology = np.array([])
        self.npoints = 0
        self.points = np.array([])
        self.ncurves = 0
        self.curves = []

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
       if not all(check) and all(check2):
           print("WARNING: Meshes not equivalent:")
           print(check)
           print(check2)

       return all(check) and all(check2)

# class definition for Geometry objects
class Geometry:

    # initialize attributes (will be filled elsewhere)
    def __init__(self,elshape=0,ndim=0,elnumnod=0,nnodes=0,nelem=0):
        self.elshape = elshape
        self.ndim = ndim
        self.elnumnod = elnumnod
        self.nnodes = nnodes
        self.nelem = nelem
        self.topology = np.array([])
        self.nodes = np.array([])

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
