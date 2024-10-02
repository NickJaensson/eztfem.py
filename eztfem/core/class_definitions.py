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


class Problem:
    """
    Define a problem with a given mesh and element degrees of freedom (DOF).

    Attributes
    ----------
    elementdof : np.ndarray
        A matrix of size (mesh.elnumnod, nvec), where nvec is the number of
        vectors of special structure.
    nphysq : int
        Number of physical quantities.
    nvec : int
        Number of vectors defined on the mesh.
    vec_elnumdegfd : np.ndarray
        Array of size (mesh.elnumnod, nvec) giving the number of DOF in each
        node of the element for the vectors.
    vec_nodnumdegfd : np.ndarray
        Array of size (mesh.nnodes+1, nvec) storing the number of DOF in each
        nodal point (accumulated) for vectors of special structure.
    vec_numdegfd : np.ndarray
        Array of length nvec storing the number of DOF for each vector of
        special structure.
    elnumdegfd : np.ndarray
        Array of size (mesh.elnumnod,) giving the number of system DOF in each
        nodal point of an element.
    nodnumdegfd : np.ndarray
        Array of length (mesh.nnodes+1) storing the number of DOF in each nodal
        point (accumulated).
    numdegfd : int
        Total number of system DOF.
    maxnoddegfd : int
        Maximum number of system DOF in the nodal points.
    maxvecnoddegfd : int
        Maximum number of vector DOF in the nodal points.

    """

    def __init__(self, mesh, elementdof, **kwargs):
        """
        Initialize the Problem instance.

        Parameters
        ----------
        mesh : object
            Mesh structure.
        elementdof : np.ndarray
            A matrix of size (mesh.elnumnod, nvec), where nvec is the number of
            vectors of special structure.
        **kwargs : dict, optional
            Optional arguments to set additional parameters.
            - nphysq : int, optional
                Number of physical quantities. Default is the number of columns
                in elementdof.

        Raises
        ------
        ValueError
            If the first dimension of elementdof does not match mesh.elnumnod.
        TypeError
            If an unexpected keyword argument is provided.

        """
        if elementdof.shape[0] != mesh.elnumnod:
            raise ValueError("First dimension of elementdof does not match \
                             mesh.elnumnod.")

        self.elementdof = elementdof

        # Optional arguments
        allowed_values = ['nphysq']
        for k in kwargs.keys():
            if k not in allowed_values:
                raise TypeError(f"Unexpected keyword argument '{k}'")

        nphysq = kwargs.get('nphysq', None)

        # Number of physical quantities
        if nphysq is None or nphysq == 0:
            self.nphysq = self.elementdof.shape[1]
        elif len(self.elementdof) < nphysq:
            raise ValueError("Length of elementdof must be at least nphysq.")
        else:
            self.nphysq = nphysq

        # Number of vectors of special structure
        self.nvec = self.elementdof.shape[1]

        # Number of degrees of freedom in each node of the element for each
        # vector of special structure
        self.vec_elnumdegfd = self.elementdof

        # Fill vec_nodnumdegfd
        self.vec_nodnumdegfd = np.zeros([mesh.nnodes+1, self.nvec], dtype=int)
        self.vec_nodnumdegfd[0, :] = 0
        self.vec_nodnumdegfd[1:, :] = -1  # Set to -1 to detect first node

        for vec in range(self.nvec):
            for elem in range(mesh.nelem):
                for node in range(mesh.elnumnod):
                    nodenr = mesh.topology[node, elem]
                    ndegfd_prev = self.vec_nodnumdegfd[nodenr + 1, vec]
                    ndegfd = self.vec_elnumdegfd[node, vec]

                    if ndegfd_prev == -1:  # First time this node
                        self.vec_nodnumdegfd[nodenr + 1, vec] = ndegfd
                    elif ndegfd != ndegfd_prev:  # Number of DOF different
                        raise ValueError(f"Different number of degrees of \
                                         freedom at nodal point {nodenr}")

        # Set isolated or inactive nodes to zero degrees of freedom
        self.vec_nodnumdegfd[self.vec_nodnumdegfd == -1] = 0

        # Accumulate vec_nodnumdegfd
        for vec in range(self.nvec):
            for nodenr in range(mesh.nnodes):
                self.vec_nodnumdegfd[nodenr + 1, vec] \
                    += self.vec_nodnumdegfd[nodenr, vec]

        # Number of degrees of freedom
        self.vec_numdegfd = self.vec_nodnumdegfd[mesh.nnodes, :]

        # Number of degrees of freedom in each node of the element
        self.elnumdegfd = np.sum(self.vec_elnumdegfd[:, :self.nphysq], axis=1)
        self.nodnumdegfd = np.sum(self.vec_nodnumdegfd[:, :self.nphysq],
                                  axis=1)
        self.numdegfd = np.sum(self.vec_numdegfd[:self.nphysq])
        self.maxnoddegfd = np.max(self.elnumdegfd)
        self.maxvecnoddegfd = np.max(self.vec_elnumdegfd)

    def __eq__(self, other):
        """
        Check equivalence with another Problem instance.

        Parameters
        ----------
        other : Problem
            Another Problem instance to compare with.

        Returns
        -------
        bool
            True if equivalent, False otherwise.

        Notes
        -----
        NOTE: see NOTE_ON_COMPARING_ARRAYS.md for the use of np.squeeze

        """
        check = [
            self.nphysq == other.nphysq,
            self.nvec == other.nvec,
            (np.squeeze(self.vec_elnumdegfd) == other.vec_elnumdegfd).all(),
            (np.squeeze(self.vec_nodnumdegfd) == other.vec_nodnumdegfd).all(),
            (self.vec_numdegfd == other.vec_numdegfd).all(),
            (self.elnumdegfd == other.elnumdegfd).all(),
            (self.nodnumdegfd == other.nodnumdegfd).all(),
            self.numdegfd == other.numdegfd,
            self.maxnoddegfd == other.maxnoddegfd,
            self.maxvecnoddegfd == other.maxvecnoddegfd
        ]

        if not all(check):
            print("WARNING: Problems not equivalent:")
            print(check)

        return all(check)


class Vector:
    """
    A class used to represent a Vector.

    Attributes
    ----------
    vec : int
        The vector value.
    u : np.ndarray
        An array to store vector data.

    Methods
    -------
    __init__(self, vec=0):
        Initializes the Vector object with the given attributes.
    __eq__(self, other):
        Checks equivalence of two Vector objects (overloads == sign).

    """

    def __init__(self, vec=0):
        """
        Initializes the Vector object with the given attributes.

        Parameters
        ----------
        vec : int or float, optional
            The vector value (default is 0).

        """
        self.vec = vec
        self.u = np.array([])

    def __eq__(self, other):
        """
        Checks equivalence of two Vector objects (overloads == sign).

        Parameters
        ----------
        other : Vector
            The other Vector object to compare with.

        Returns
        -------
        bool
            True if the Vector objects are equivalent, False otherwise.

        Notes
        -----
        See NOTE_ON_COMPARING_ARRAYS.md for the use of np.squeeze.

        """
        check = [self.vec == other.vec,
                 np.allclose(np.squeeze(self.u), other.u, atol=1e-12, rtol=0)]
        return all(check)


class User:
    """
    Define the User class: this class is mainly a collection of parameters and
    data that should be accessible to the element routines.

    Attributes
    ----------
    wg : np.ndarray
        Weights of Gauss points.
    xr : np.ndarray
        Location of Gauss points.
    phi : np.ndarray
        Values of phi shape function in Gauss points.
    dphi : np.ndarray
        Values of phi shape function derivatives in Gauss points.
    psi : np.ndarray
        Values of psi shape function in Gauss points.

    Methods
    -------
    __init__(self):
        Initializes the User object with the given attributes.
    __eq__(self, other):
        Checks equivalence of two User objects (overloads == sign).

    """

    def __init__(self):
        """
        Initializes the User object with empty numpdy.ndarray.
        """
        self.wg = np.array([])
        self.xr = np.array([])
        self.phi = np.array([])
        self.dphi = np.array([])
        self.psi = np.array([])

    def __eq__(self, other):
        """
        Checks equivalence of two User objects (overloads == sign).

        Parameters
        ----------
        other : User
            The other User object to compare with.

        Returns
        -------
        bool
            True if the User objects are equivalent, False otherwise.

        Notes
        -----
        See NOTE_ON_COMPARING_ARRAYS.md for the use of np.squeeze.
        Only a selection of possible attributes is checked.

        """
        # Get list of names of attributes
        attributes_self = [attr for attr in self.__dict__.keys()]
        attributes_other = [attr for attr in other.__dict__.keys()]

        # Check if attributes exist in both Users
        check0 = [
            ('wg' in attributes_self) is ('wg' in attributes_other),
            ('xr' in attributes_self) is ('xr' in attributes_other),
            ('phi' in attributes_self) is ('phi' in attributes_other),
            ('psi' in attributes_self) is ('psi' in attributes_other),
            ('dphi' in attributes_self) is ('dphi' in attributes_other),
            ('coorsys' in attributes_self) is ('coorsys' in attributes_other),
            ('alpha' in attributes_self) is ('alpha' in attributes_other),
            ('mu' in attributes_self) is ('mu' in attributes_other),
            ('funcnr' in attributes_self) is ('funcnr' in attributes_other)
        ]

        if not all(check0):
            print("WARNING: Users do not have equal attributes:")
            print(check0)
            return False

        # Check if existing attributes have same values
        check1 = [True]  # Avoid empty list
        if 'wg' in attributes_self:
            check1.append(np.allclose(self.wg, other.wg, atol=1e-15, rtol=0))
        if 'xr' in attributes_self:
            check1.append(np.allclose(np.squeeze(self.xr), other.xr,
                                      atol=1e-15, rtol=0))
        if 'phi' in attributes_self:
            check1.append(np.allclose(np.squeeze(self.phi), other.phi,
                                      atol=1e-15, rtol=0))
        if 'psi' in attributes_self:
            check1.append(np.allclose(np.squeeze(self.psi), other.psi,
                                      atol=1e-15, rtol=0))
        if 'dphi' in attributes_self:
            check1.append(np.allclose(self.dphi, other.dphi, atol=1e-15,
                                      rtol=0))
        if 'coorsys' in attributes_self:
            check1.append(self.coorsys == other.coorsys)
        if 'alpha' in attributes_self:
            check1.append(self.alpha == other.alpha)
        if 'mu' in attributes_self:
            check1.append(self.mu == other.mu)
        if 'funcnr' in attributes_self:
            check1.append(self.funcnr == other.funcnr)

        # Print a warning when not equivalent (for debugging purposes)
        if not all(check1):
            print("WARNING: User attributes not equivalent:")
            print(check1)

        return all(check1)
