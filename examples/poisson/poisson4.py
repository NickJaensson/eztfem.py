# Poisson problem on a unit square with Dirichlet and natural boundary
# conditions. Manufactured solution.

import numpy as np
import eztfem as ezt
from func import func
from scipy.sparse.linalg import spsolve


def main():

    # create mesh

    print('mesh')
    mesh = ezt.quadrilateral2d([20, 20], 'quad9', length=np.array([1.2, 1]))

    # define the problem

    print('problem_definition')
    elementdof = np.array([[1, 1, 1, 1, 1, 1, 1, 1, 1],
                           [2, 2, 2, 2, 2, 2, 2, 2, 2]], dtype=int).transpose()
    problem = ezt.Problem(mesh, elementdof, nphysq=1)

    # define Gauss integration and basis functions

    user = ezt.User()
    shape = 'quad'

    print('gauss_legendre')
    user.xr, user.wg = ezt.gauss_legendre(shape, n=3)

    print('basis_function phi')
    user.phi, user.dphi = ezt.basis_function(shape, 'Q2', user.xr)

    # user struct for setting problem coefficients, ...

    user.coorsys = 0
    user.alpha = 1
    user.funcnr = 7
    user.func = func

    # assemble the system matrix and vector

    print('build_system')
    A, f = ezt.build_system(mesh, problem, ezt.poisson_elem, user)

    # define Gauss integration and basis functions (for boundary integral)

    print('gauss_legendre')
    [xr, user.wg] = ezt.gauss_legendre('line', n=3)

    print('basis_function phi')
    [user.phi, user.dphi] = ezt.basis_function('line', 'P2', xr)

    # add natural boundary condition

    print('add_boundary_elements')
    user.funcnr = 8
    ezt.add_boundary_elements(mesh, problem, f, ezt.poisson_natboun_curve,
                              user, curve=1)

    # define essential boundary conditions (Dirichlet)

    print('define_essential')
    iess = ezt.define_essential(mesh, problem, 'curves', [0, 2, 3])

    # fill values for the essential boundary conditions

    print('fill_system_vector')
    uess = ezt.fill_system_vector(mesh, problem, 'curves', [0, 2, 3],
                                  func, funcnr=6)

    # apply essential boundary conditions to the system

    print('apply_essential')
    ezt.apply_essential(A, f, uess, iess)

    # solve the system

    print('solve')
    u = spsolve(A.tocsr(), f)

    # compare with exact solution

    print('Difference with exact solution:')
    uex = ezt.fill_system_vector(mesh, problem, 'nodes',
                                 np.arange(mesh.nnodes), func, funcnr=6)

    maxdiff = max(abs(u-uex))
    print(maxdiff)

    # gradient (dudx,dudy) of the solution

    print('gradient (dudx,dudy)')
    xr = ezt.refcoor_nodal_points(mesh)
    [user.phi, user.dphi] = ezt.basis_function('quad', 'Q2', xr)
    user.u = u

    gradu = ezt.deriv_vector(mesh, problem, ezt.poisson_deriv, user)

    return maxdiff


if __name__ == '__main__':
    main()
