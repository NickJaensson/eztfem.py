import numpy as np

def pos_array(problem, nodes, **kwargs):
    """
    Get the index of the system degrees of freedom in the given nodes.

    Parameters:
    problem: object representing the problem structure
    nodes: An array of node numbers
    **kwargs: Optional arguments
        - physq: Array of physical quantity numbers, default: all physical quantities
        - order: The sequence order of the degrees of freedom in the nodes
    
    Returns:
    pos: List of arrays containing the positions of the degrees of freedom of each physq
    ndof: List of the number of degrees in pos of each physq
    """

    # Set default optional arguments
    physq = np.arange(problem.nphysq)
    order = 'DN'
    
    # Override optional arguments if provided
    if 'physq' in kwargs:
        physq = kwargs['physq']
    if 'order' in kwargs:
        order = kwargs['order']

    # Convert physq to a list if an int is supplied
    if isinstance(physq, (int, np.integer)):
        physq = [physq]

    # Convert nodes to a numpy array if an int is supplied
    if isinstance(nodes, (int, np.integer)):
        nodes = np.array([nodes])

    pos = [None] * len(physq)
    ndof = np.zeros(len(physq))
    lpos = np.zeros(problem.maxnoddegfd * nodes.shape[0],dtype=int)
    if order == 'ND':
        for i, phq in enumerate(physq):
            dof = 0
            for nodenr in nodes:
                bp = problem.nodnumdegfd[nodenr] + \
                     sum(problem.vec_nodnumdegfd[nodenr+1, :phq] - problem.vec_nodnumdegfd[nodenr, :phq])
                nndof = problem.vec_nodnumdegfd[nodenr+1, phq] - problem.vec_nodnumdegfd[nodenr, phq]
                lpos[dof:dof+nndof] = np.arange(bp, bp+nndof)
                dof += nndof
            pos[i] = lpos[:dof].tolist() # convert to list: apparently mixing np arrays and lists goes wrong!
            ndof[i] = dof
    elif order == 'DN':
        for i, phq in enumerate(physq):
            maxdeg = max(problem.vec_nodnumdegfd[nodes+1, phq] - problem.vec_nodnumdegfd[nodes, phq])
            dof = 0
            for deg in range(maxdeg):
                for nodenr in nodes:
                    bp = problem.nodnumdegfd[nodenr] + \
                         sum(problem.vec_nodnumdegfd[nodenr+1, :phq] - problem.vec_nodnumdegfd[nodenr, :phq])
                    nndof = problem.vec_nodnumdegfd[nodenr+1, phq] - problem.vec_nodnumdegfd[nodenr, phq]
                    if deg < nndof:
                        lpos[dof] = bp + deg
                        dof += 1
            pos[i] = lpos[:dof].tolist() # convert to list: apparently mixing np arrays and lists goes wrong!
            ndof[i] = dof

    return pos, ndof
