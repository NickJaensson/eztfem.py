import numpy as np

class Problem:
    
#PROBLEM_DEFINITION  Define the problem
#  problem = PROBLEM_DEFINITION ( mesh, elementdof, 'optarg1', value1, ... )
#  input:
#    mesh: mesh structure
#    elementdof: a matrix of size mesh.elnumnod x nvec, where nvec
#                is the number of "vectors of special structure".
#                Each column of the matrix gives the number of degrees of
#                freedom in each node of the element for the vectors.
#  optional arguments:
#    string, value couples to set optional parameters:
#    'nphysq', the number of "physical quantities". The first nphysq vectors
#              of special structure (first nphysq columns of elementdof) define
#              the physical quantities.
#              default: size(elementdof,1)
#  output:
#    problem: the problem structure, having the components
#       nphysq: number of physical quantities
#       nvec: number of vectors defined on the mesh. The physical quantities
#             are defined by the first nphysq vectors.
#       vec_elnumdegfd: an array of size (mesh.elnumnod,nvec), where each
#             column of the matrix gives the number of degrees of freedom in
#             each node of the element for the vectors.
#       vec_nodnumdegfd: an array of length mesh.nnodes+1 x nvec storing the
#             number of degrees of freedom in each nodal point (accumulated)
#             of vectors of special structure as follows:
#                  vec_nodnumdegfd(1,vec) = 0  
#             number of degrees of freedom for nodal point n is
#                  vec_nodnumdegfd(n+1,vec) - vec_nodnumdegfd(n,vec)
#             for vector of special structure vec
#       vec_numdegfd: an array of length nvec storing the number of degrees
#             of freedom for each vector of special structure
#       elnumdegfd: an array of size mesh.elnumnod, giving the number of
#             (system) degrees of freedom in each nodal point of an element
#       nodnumdegfd: an array of length mesh.nnodes+1 storing the
#             number of degrees of freedom in each nodal point (accumulated)
#                  nodnumdegfd(1) = 0  
#             number of degrees of freedom for nodal point n is
#                  nodnumdegfd(n+1) - vec_nodnumdegfd(n)
#       numdegfd: the number of (system) degrees of freedom 
#       maxnoddegfd: the maximum number of (system) degrees of freedom in the
#                    the nodal points 
#       maxvecnoddegfd: the maximum number of vector degrees of freedom in the
#                    the nodal points 

    def __init__(self, mesh, elementdof, **kwargs):

        # Check if the first dimension of elementdof matches mesh.elnumnod
        if len(elementdof[0]) != mesh.elnumnod:
            raise ValueError("First dimension of elementdof does not match elnumnod.")
        
        self.elementdof = elementdof

        # Optional arguments
        nphysq = kwargs.get('nphysq', None)

        # Number of physical quantities
        if not nphysq:
            self.nphysq = self.elementdof.shape[1]
        elif len(self.elementdof) < nphysq:
            raise ValueError("Length of elementdof must be at least nphysq.")
        else:
            self.nphysq = nphysq
            
        # Number of vectors of special structure
        self.nvec = len(self.elementdof)

        # Number of degrees of freedom in each node of the element for each vector
        self.vec_elnumdegfd = self.elementdof

        # fill nodnumdegfd
        self.vec_nodnumdegfd = np.zeros([mesh.nnodes+1,self.nvec],dtype=int)
        self.vec_nodnumdegfd[0,:]  = 0
        self.vec_nodnumdegfd[1:,:] = -1  # set to -1 to detect first node

        for vec in range(self.nvec):
            for elem in range(mesh.nelem):
                for node in range(mesh.elnumnod):
                    nodenr = mesh.topology[node, elem]
                    ndegfd_prev = self.vec_nodnumdegfd[nodenr + 1, vec]
                    ndegfd = self.vec_elnumdegfd[vec,node]

                    if ndegfd_prev == -1:
                        self.vec_nodnumdegfd[nodenr + 1, vec] = ndegfd
                    elif ndegfd != ndegfd_prev:
                        raise ValueError(f"Different number of degrees of freedom at nodal point {nodenr}")

        # set isolated or inactive nodes to zero degrees of freedom
        self.vec_nodnumdegfd[self.vec_nodnumdegfd == -1] = 0

        # accumulate vec.nodnumdegfd
        for vec in range(self.nvec):
            for nodenr in range(mesh.nnodes):
                self.vec_nodnumdegfd[nodenr + 1, vec] += self.vec_nodnumdegfd[nodenr, vec]


        # number of degrees of freedom
        self.vec_numdegfd = self.vec_nodnumdegfd[mesh.nnodes,:]


        # Now all degrees in the nodes for the system vector        

        # number of degrees of freedom in each node of the element
        self.elnumdegfd = np.sum(self.vec_elnumdegfd[:, :self.nphysq], axis=1)
        self.nodnumdegfd = np.sum(self.vec_nodnumdegfd[:, :self.nphysq], axis=1)
        self.numdegfd = np.sum(self.vec_numdegfd[:self.nphysq])
        self.maxnoddegfd = np.max(self.elnumdegfd)
        self.maxvecnoddegfd = np.max(self.vec_elnumdegfd)


