'''
Poisson problem on a unit square with Dirichlet and natural boundary
conditions.
Manufactured solution.
'''

import numpy as np
from scipy.sparse.linalg import spsolve
from func import func
import eztfem as ezt


def main():
    """
    Returns
    -------
    max_diff : float
        Maximum difference of numerical and analytical solution.
    mesh : Mesh
        Finite element mesh.
    problem : Problem
        Stokes problem definition.
    u : numpy.ndarray
        Nodal values of u (solution vector).
    derivs : dict
        Dictionary of derived fields (Vector instances).
    """

    # create mesh

    mesh = ezt.quadrilateral2d([20, 20], 'quad9', length=np.array([1.2, 1]))

    # define the problem

    elementdof = np.array([[1, 1, 1, 1, 1, 1, 1, 1, 1],
                           [2, 2, 2, 2, 2, 2, 2, 2, 2]], dtype=int).transpose()
    problem = ezt.Problem(mesh, elementdof, nphysq=1)

    # define Gauss integration and basis functions

    user = ezt.User()
    shape = 'quad'

    user.xr, user.wg = ezt.gauss_legendre(shape, n=3)
    user.phi, user.dphi = ezt.basis_function(shape, 'Q2', user.xr)

    # user struct for setting problem coefficients, ...

    user.coorsys = 0
    user.alpha = 1
    user.funcnr = 7
    user.func = func

    # assemble the system matrix and vector

    system_matrix, rhs = ezt.build_system(mesh, problem, ezt.poisson_elem,
                                          user)

    # define Gauss integration and basis functions (for boundary integral)

    [xr, user.wg] = ezt.gauss_legendre('line', n=3)
    [user.phi, user.dphi] = ezt.basis_function('line', 'P2', xr)

    # add natural boundary condition

    user.funcnr = 8
    ezt.add_boundary_elements(mesh, problem, rhs, ezt.poisson_natboun_curve,
                              user, curve=1)

    # define essential boundary conditions (Dirichlet)

    iess = ezt.define_essential(mesh, problem, 'curves', [0, 2, 3])

    # fill values for the essential boundary conditions

    uess = ezt.fill_system_vector(mesh, problem, 'curves', [0, 2, 3],
                                  func, funcnr=6)

    # apply essential boundary conditions to the system

    ezt.apply_essential(system_matrix, rhs, uess, iess)

    # solve the system

    u = spsolve(system_matrix.tocsr(), rhs)

    # compare with exact solution

    uex = ezt.fill_system_vector(mesh, problem, 'nodes',
                                 np.arange(mesh.nnodes), func, funcnr=6)

    max_diff = max(abs(u-uex))

    # gradient (dudx,dudy) of the solution

    xr = ezt.refcoor_nodal_points(mesh)
    [user.phi, user.dphi] = ezt.basis_function('quad', 'Q2', xr)
    user.u = u

    derivs = {}
    derivs["grad_u"] = ezt.deriv_vector(mesh, problem, ezt.poisson_deriv, user)

    return max_diff, mesh, problem, u, derivs


if __name__ == '__main__':
    results = main()
    print('max diff = ', results[0])
