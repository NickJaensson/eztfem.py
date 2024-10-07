import numpy as np
from .pos_array import pos_array


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

    def __init__(self, mesh, elementdof, nphysq=None):
        """
        Initializes the Problem object for a given mesh, elementdof and
        nphys (optional). All other attributes are filled in the initialization
        function.

        Parameters for initialization
        -----------------------------
        mesh : Mesh
            Mesh object.
        elementdof : np.ndarray
            A matrix of size (mesh.elnumnod, nvec), where nvec is the number of
            vectors of special structure.
        nphysq : int, optional
            Number of physical quantities. Default is the number of columns
            in elementdof.

        """
        if elementdof.shape[0] != mesh.elnumnod:
            raise ValueError("First dimension of elementdof does not match \
                             mesh.elnumnod.")

        self.elementdof = elementdof

        # Number of physical quantities
        if nphysq is None:
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


def define_essential(mesh, problem, geometry, numbers, **kwargs):
    """
    Get the indices of the essential degrees of freedom.

    Parameters
    ----------
    mesh : Mesh
        Mesh object.
    problem : Problem
        Problem object.
    geometry : str
        The 'geometry' to put essential degrees on:
        'nodes' : in the nodes given by numbers.
        'points' : in the points given by numbers.
        'curves' : in the curves given by numbers.
    numbers : array_like
        An array of 'geometry' numbers.

    Keyword arguments
    -----------------
    physq : int, optional
        Physical quantity number. Default is 0.
    degfd : int, optional
        Degree of freedom within the physical quantity. Default is 0.
    iessp : int, optional
        Index of previously defined essential degrees. Newly defined
        boundary conditions will be added. Default is 0.

    Returns
    -------
    iess : numpy.ndarray
        Index of defined essential degrees.

    """

    # Optional arguments
    physq = kwargs.get('physq', 0)
    degfd = kwargs.get('degfd', 0)
    iessp = kwargs.get('iessp', 0)
    add = 'iessp' in kwargs

    # Convert numbers to a list if an int is supplied
    if isinstance(numbers, (int, np.integer)):
        numbers = [numbers]

    # Define geometry
    if geometry == 'nodes':
        nodes = numbers
        nnodes = len(numbers)
    elif geometry == 'points':
        nodes = mesh.points[numbers]
        nnodes = len(numbers)
    elif geometry == 'curves':
        nnodes = sum([mesh.curves[curve].nnodes for curve in numbers])
        nodes = np.array([], dtype=int)
        for curve in numbers:
            nodes = np.append(nodes, mesh.curves[curve].nodes)
    else:
        raise ValueError(f"Invalid geometry: {geometry}")

    # Determine pos array
    pos = np.zeros(nnodes, dtype=int)
    ipos = 0

    for node in nodes:
        posn, ndof = pos_array(problem, node, physq=physq, order='ND')

        if degfd > ndof[0]:
            continue

        pos[ipos] = posn[0][degfd]  # Python 0-based indexing
        ipos += 1

    # Output iess
    if add:
        # Convert iessp to a numpy array if an int is supplied
        if isinstance(iessp, (int, np.integer)):
            iessp = np.array([iessp])
        iess = np.unique(np.concatenate((pos[:ipos], iessp)))
    else:
        iess = np.unique(pos[:ipos])

    return iess
