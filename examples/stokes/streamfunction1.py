'''
Compute streamfunction from a velocity field.
'''

import numpy as np
from scipy.sparse.linalg import spsolve
import eztfem as ezt


def main(mesh=None, problem=None, u=None):
    """
    Returns
    -------
    max_streamfunction : float
        Maximum streamfunction value.
    streamfunction : numpy.ndarray
        Nodal values of the streamfunction.
    """

    if mesh is None:
        raise ValueError("Error streamfunction1: mesh should be present")

    if problem is None:
        raise ValueError("Error streamfunction1: problem should be present")

    if u is None:
        raise ValueError("Error streamfunction1: u should be present")

    # define the problem

    elementdof = np.array([[1, 1, 1, 1, 1, 1, 1, 1, 1],
                           [2, 2, 2, 2, 2, 2, 2, 2, 2]], dtype=int).transpose()
    problem_s = ezt.Problem(mesh, elementdof, nphysq=1)

    # define Gauss integration and basis functions

    user = ezt.User()
    shape = 'quad'
    user.xr, user.wg = ezt.gauss_legendre(shape, n=3)
    user.phi, user.dphi = ezt.basis_function(shape, 'Q2', user.xr)

    # user struct for setting problem coefficients, ...

    user.coorsys = 0

    # Get velocities node for node
    pos, _ = ezt.pos_array(problem, np.arange(mesh.nnodes), physq=0,
                           order='ND')
    user.v = u[pos[0]]

    # assemble the system matrix and vector

    system_matrix, rhs = ezt.build_system(mesh, problem_s,
                                          ezt.streamfunction_elem, user,
                                          posvectors=True)

    # define Gauss integration and basis functions (for boundary integral)

    [xr, user.wg] = ezt.gauss_legendre('line', n=3)
    [user.phi, user.dphi] = ezt.basis_function('line', 'P2', xr)

    # add natural boundary condition

    for crv in range(3):
        ezt.add_boundary_elements(mesh, problem_s, rhs,
                                  ezt.streamfunction_natboun_curve, user,
                                  posvectors=True, curve=crv)

    # define essential boundary conditions (Dirichlet)

    iess = ezt.define_essential(mesh, problem_s, 'points', [0])

    # fill values for the essential boundary conditions

    uess = np.zeros([problem_s.numdegfd])

    # apply essential boundary conditions to the system

    ezt.apply_essential(system_matrix, rhs, uess, iess)

    # solve the system

    streamfunction = spsolve(system_matrix.tocsr(), rhs)

    return max(streamfunction), streamfunction


if __name__ == '__main__':

    import stokes1

    _, mesh_stokes1, problem_stokes1, u_stokes1, _ = stokes1.main()

    results = main(mesh_stokes1, problem_stokes1, u_stokes1)

    print('Maximum value ', results[0])
