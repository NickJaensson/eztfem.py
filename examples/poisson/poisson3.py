'''
Poisson problem on an L-shape domain with Dirichlet boundary conditions
This problem solves Example 1.1.2 of the book of Elman et al.
'''

import numpy as np
from scipy.sparse.linalg import spsolve
from func import func
import eztfem as ezt


def main():
    """
    Returns
    -------
    max_u : float
        Maximum solution value.
    mesh : Mesh
        Finite element mesh.
    problem : Problem
        Stokes problem definition.
    u : numpy.ndarray
        Nodal values of solution u.
    """

    # create mesh

    mesh = ezt.l_shape2d([20, 20, 20, 20], 'quad4',
                         factor=np.array([10, 10, 10, 10]))

    # define the problem

    elementdof = np.array([[1, 1, 1, 1],
                           [2, 2, 2, 2]], dtype=int).transpose()
    problem = ezt.Problem(mesh, elementdof, nphysq=1)

    # define Gauss integration and basis functions

    user = ezt.User()
    shape = 'quad'

    user.xr, user.wg = ezt.gauss_legendre(shape, n=2)
    user.phi, user.dphi = ezt.basis_function(shape, 'Q1', user.xr)

    # user struct for setting problem coefficients, ...

    user.coorsys = 0
    user.alpha = 1
    user.funcnr = 5
    user.func = func

    # assemble the system matrix and vector

    system_matrix, rhs = ezt.build_system(mesh, problem, ezt.poisson_elem,
                                          user)

    # define essential boundary conditions (Dirichlet)

    iess = ezt.define_essential(mesh, problem, 'curves',
                                [0, 1, 2, 3, 4, 5, 6, 7])

    # fill values for the essential boundary conditions

    uess = np.zeros([problem.numdegfd])

    # apply essential boundary conditions to the system

    ezt.apply_essential(system_matrix, rhs, uess, iess)

    # solve the system

    u = spsolve(system_matrix.tocsr(), rhs)

    return max(u), mesh, problem, u


if __name__ == '__main__':
    results = main()
    print('max u = ', results[0])
