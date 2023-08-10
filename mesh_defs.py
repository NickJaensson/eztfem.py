import numpy as np

class Mesh:
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

class Geometry:
    def __init__(self,elshape=0,ndim=0,elnumnod=0,nnodes=0,nelem=0):
        self.elshape = elshape
        self.ndim = ndim
        self.elnumnod = elnumnod
        self.nnodes = nnodes
        self.nelem = nelem
        self.topology = np.array([])
        self.nodes = np.array([])