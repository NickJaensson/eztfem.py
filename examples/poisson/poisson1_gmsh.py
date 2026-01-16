'''
Poisson problem on a unit square with Dirichlet boundary conditions
Manufactured solution.
Mesh generated with Gmsh (triangular elements of order 2).
'''

import numpy as np
from scipy.sparse.linalg import spsolve
from func import func
import gmsh
import eztfem as ezt


import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))

from gmsh_meshes import unit_square


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

    # set global Gmsh options

    gmsh.initialize()
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.08)
    gmsh.option.setNumber("Mesh.ElementOrder", 2)
    gmsh.option.setNumber("Mesh.HighOrderOptimize", 1)

    # create mesh

    mesh = ezt.gmsh_mesh2d(model_builder=unit_square.build_unit_square)

    # define the problem

    elementdof = np.array([[1, 1, 1, 1, 1, 1],
                           [2, 2, 2, 2, 2, 2]], dtype=int).transpose()
    problem = ezt.Problem(mesh, elementdof, nphysq=1)

    # define Gauss integration and basis functions

    user = ezt.User()
    shape = 'triangle'

    user.xr, user.wg = ezt.gauss_legendre(shape, p=3)
    user.phi, user.dphi = ezt.basis_function(shape, 'P2', user.xr)

    # user struct for setting problem coefficients, ...

    user.coorsys = 0
    user.alpha = 1
    user.funcnr = 4
    user.func = func

    # assemble the system matrix and vector

    system_matrix, rhs = ezt.build_system(mesh, problem, ezt.poisson_elem,
                                          user)

    # define essential boundary conditions (Dirichlet)

    iess = ezt.define_essential(mesh, problem, 'curves', [0, 1, 2, 3])

    # fill values for the essential boundary conditions

    uess = ezt.fill_system_vector(mesh, problem, 'curves', [0, 1, 2, 3], func,
                                  funcnr=3)

    # apply essential boundary conditions to the system

    ezt.apply_essential(system_matrix, rhs, uess, iess)

    # solve the system

    u = spsolve(system_matrix.tocsr(), rhs)

    # compare with exact solution

    uex = ezt.fill_system_vector(mesh, problem, 'nodes',
                                 np.arange(mesh.nnodes), func, funcnr=3)

    max_diff = max(abs(u-uex))

    # gradient (dudx,dudy) of the solution

    xr = ezt.refcoor_nodal_points(mesh)
    [user.phi, user.dphi] = ezt.basis_function('triangle', 'P2', xr)
    user.u = u

    derivs = {}
    derivs["grad_u"] = ezt.deriv_vector(mesh, problem, ezt.poisson_deriv, user)

    return max_diff, mesh, problem, u, derivs


if __name__ == '__main__':
    results = main()
    print('max diff = ', results[0])
