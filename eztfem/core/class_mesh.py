import numpy as np


class Mesh:
    """
    A class used to represent a Mesh.

    Attributes
    ----------
    ndim : int
       Dimension of space (ndim=1 or 2).
    nnodes : int
        Number of nodes.
    nelem : int
        Number of elements.
    elshape : int
        Shape number of the elements. Table of shapes:
        - 1 : 2-node line elements.
        - 2 : 3-node line elements.
        - 3 : 3-node triangle.
        - 4 : 6-node triangle.
        - 5 : 4-node quadrilateral.
        - 6 : 9-node quadrilateral.
        - 7 : 7-node triangle.
        - 9 : 5-node quadrilateral.
        - 10 : 4-node triangle.
    elnumnod : int
        Number of nodes in a single element.
    npoints : int
        Number of points.
    ncurves : int
        Number of curves.
    topology : np.ndarray
        Array of size (elnumnod, nelem), where topology[:, elem] is an
        array of global node numbers element elem is connected to.
    coor : np.ndarray
        Array of size (nnodes, ndim), where coor[i, :] are the coordinates
        of node i, with 1 <= i <= nnodes.
    points : np.ndarray
        An array containing the node numbers of the points, hence points[i]
        is the node of point i, with 1 <= i <= npoints.
    curves : list
        An array of objects of type Geometry

    Methods
    -------
    __init__(self, ndim=0, nnodes=0, nelem=0, elshape=0, elnumnod=0,
        npoints=0, ncurves=0, topology=None, coor=None, points=None,
        curves=None):
        Initializes the Mesh object with the given attributes.
    __eq__(self, other):
        Checks equivalence of two Mesh objects (overloads == sign)

    """

    def __init__(self, ndim=0, nnodes=0, nelem=0, elshape=0, elnumnod=0,
                 npoints=0, ncurves=0, topology=None, coor=None, points=None,
                 curves=None):
        """
        Initializes the Mesh object with the given attributes.

        Parameters
        ----------
        ndim : int, optional
            Number of dimensions (default is 0).
        nnodes : int, optional
            Number of nodes (default is 0).
        nelem : int, optional
            Number of elements (default is 0).
        elshape : int, optional
            Shape of the elements (default is 0).
        elnumnod : int, optional
            Number of nodes per element (default is 0).
        npoints : int, optional
            Number of points (default is 0).
        ncurves : int, optional
            Number of curves (default is 0).
        topology : np.ndarray, optional
            Topology of the mesh (default is None).
        coor : np.ndarray, optional
            Coordinates of the nodes (default is None).
        points : np.ndarray, optional
            Points in the mesh (default is None).
        curves : list, optional
            Curves in the mesh (default is None).

        Notes
        -----
        It is not recommended to use mutable objects as default values in
        function def (e.g., for coor, points or curves). See:
        docs.python.org/3/tutorial/controlflow.html#default-argument-values

        """
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

    def __eq__(self, other):
        """
        Checks equivalence of two Mesh objects (overloads == sign).

        Parameters
        ----------
        other : Mesh
            The other Mesh object to compare with.

        Returns
        -------
        bool
            True if the Mesh objects are equivalent, False otherwise.

        Notes
        -----
        NOTE: see NOTE_ON_COMPARING_ARRAYS.md for the use of np.squeeze

        """
        check = [self.ndim == other.ndim,
                 self.nnodes == other.nnodes,
                 np.allclose(np.squeeze(self.coor), other.coor,
                             atol=1e-15, rtol=0.0),
                 self.nelem == other.nelem,
                 self.elshape == other.elshape,
                 self.elnumnod == other.elnumnod,
                 (np.squeeze(self.topology) == other.topology).all(),
                 self.npoints == other.npoints,
                 (self.points == other.points).all(),
                 self.ncurves == other.ncurves]
        check2 = [self.curves[i] == other.curves[i]
                  for i in range(self.ncurves)]

        if not (all(check) and all(check2)):
            print("WARNING: Meshes not equivalent:")
            print(check)  # print check for debugging purposes
            print(check2)  # print check for debugging purposes

        return all(check) and all(check2)


class Geometry:
    """
    A class used to represent a Geometry.

    Attributes
    ----------
    elshape : int
        Shape number of the elements. Table of shapes:
        - 1 : 2-node elements.
        - 2 : 3-node elements.
    ndim : int
        Dimension of space (ndim=2).
    elnumnod : int
        Number of nodes in a single element.
    nnodes : int
        Number of nodes.
    nelem : int
        Number of elements.
    topology : np.ndarray
        Array of size (elnumnod, nelem, 2), where topology[:, elem, 0]
        is an array of local node numbers element elem is connected to,
        and topology[:, elem, 1] is an array of global node numbers
        element elem is connected to.
    nodes : np.ndarray
        Array of size nnodes containing the global node numbers, i.e.,
        nodes[i] is the global node number of local node i.

    Methods
    -------
    __init__(self, elshape=0, ndim=0, elnumnod=0, nnodes=0, nelem=0,
        topology=None, nodes=None):
        Initializes the Geometry object with the given attributes.
    __eq__(self, other):
        Checks equivalence of two Geometry objects (overloads == sign).

    """

    def __init__(self, elshape=0, ndim=0, elnumnod=0, nnodes=0, nelem=0,
                 topology=None, nodes=None):
        """
        Initializes the Geometry object with the given attributes.

        Parameters
        ----------
        elshape : int, optional
            Shape of the elements (default is 0).
        ndim : int, optional
            Number of dimensions (default is 0).
        elnumnod : int, optional
            Number of nodes per element (default is 0).
        nnodes : int, optional
            Number of nodes (default is 0).
        nelem : int, optional
            Number of elements (default is 0).
        topology : np.ndarray, optional
            Topology of the geometry (default is None).
        nodes : np.ndarray, optional
            Nodes in the geometry (default is None).

        """
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

    def __eq__(self, other):
        """
        Checks equivalence of two Geometry objects (overloads == sign).

        Parameters
        ----------
        other : Geometry
            The other Geometry object to compare with.

        Returns
        -------
        bool
            True if the Geometry objects are equivalent, False otherwise.

        Notes
        -----
        NOTE: see NOTE_ON_COMPARING_ARRAYS.md for the use of np.squeeze

        """
        check = [self.ndim == other.ndim,
                 self.elshape == other.elshape,
                 self.elnumnod == other.elnumnod,
                 self.nnodes == other.nnodes,
                 self.nelem == other.nelem,
                 (self.topology == other.topology).all(),
                 (self.nodes == other.nodes).all()]

        if not all(check):
            print("WARNING: Geometries not equivalent:")
            print(check)  # print check for debugging purposes

        return all(check)
