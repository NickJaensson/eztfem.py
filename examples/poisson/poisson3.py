# Poisson problem on an L-shape domain with Dirichlet boundary conditions
# This problem solves Example 1.1.2 of the book of Elman et al.

import numpy as np
import eztfem as ezt
from func import func
from scipy.sparse.linalg import spsolve


def main():

    # create mesh

    print('mesh')
    mesh = ezt.l_shape2d([20, 20, 20, 20], 'quad4',
                         factor=np.array([10, 10, 10, 10]))

    # define the problem

    print('problem_definition')
    elementdof = np.array([[1, 1, 1, 1],
                           [2, 2, 2, 2]], dtype=int).transpose()
    problem = ezt.Problem(mesh, elementdof, nphysq=1)

    # define Gauss integration and basis functions

    user = ezt.User()
    shape = 'quad'

    print('gauss_legendre')
    user.xr, user.wg = ezt.gauss_legendre(shape, n=2)

    print('basis_function phi')
    user.phi, user.dphi = ezt.basis_function(shape, 'Q1', user.xr)

    # user struct for setting problem coefficients, ...

    user.coorsys = 0
    user.alpha = 1
    user.funcnr = 5
    user.func = func

    # assemble the system matrix and vector

    print('build_system')
    A, f = ezt.build_system(mesh, problem, ezt.poisson_elem, user)

    # define essential boundary conditions (Dirichlet)

    print('define_essential')
    iess = ezt.define_essential(mesh, problem, 'curves',
                                [0, 1, 2, 3, 4, 5, 6, 7])

    # fill values for the essential boundary conditions

    print('fill_system_vector')
    uess = np.zeros([problem.numdegfd])

    # apply essential boundary conditions to the system

    print('apply_essential')
    ezt.apply_essential(A, f, uess, iess)

    # solve the system

    print('solve')
    u = spsolve(A.tocsr(), f)

    # maximum value

    print('Maximum value ', max(u))

    return max(u)


if __name__ == '__main__':
    main()
