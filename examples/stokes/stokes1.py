'''
Stokes problem on a unit square with Dirichlet boundary conditions.
Lid-driven cavity flow.
'''

import numpy as np
from scipy.sparse.linalg import spsolve
from func import func
import eztfem as ezt


def main():
    """
    Returns
    -------
    max_omega : float
        Maximum vorticity magnitude.
    mesh : Mesh
        Finite element mesh.
    problem : Problem
        Stokes problem definition.
    u : numpy.ndarray
        Nodal values of velocity and pressure (solution vector).
    derivs : dict
        Dictionary of derived fields (Vector instances).
    """

    # create mesh

    mesh = ezt.quadrilateral2d([20, 20], 'quad9')

    # define the problem

    elementdof = np.array([[2, 2, 2, 2, 2, 2, 2, 2, 2],
                           [1, 0, 1, 0, 1, 0, 1, 0, 0],
                           [1, 1, 1, 1, 1, 1, 1, 1, 1]], dtype=int).transpose()
    problem = ezt.Problem(mesh, elementdof, nphysq=2)

    # define Gauss integration and basis functions

    user = ezt.User()
    shape = 'quad'

    user.xr, user.wg = ezt.gauss_legendre(shape, n=3)
    user.phi, user.dphi = ezt.basis_function(shape, 'Q2', user.xr)
    user.psi, _ = ezt.basis_function(shape, 'Q1', user.xr)

    # user struct for setting problem coefficients, ...

    user.coorsys = 0
    user.mu = 1
    user.funcnr = 0
    user.func = func  # not used when funcnr == 0

    # assemble the system matrix and right hand side

    system_matrix, rhs = ezt.build_system(mesh, problem, ezt.stokes_elem, user)

    # define essential boundary conditions (Dirichlet)

    iess = ezt.define_essential(mesh, problem, 'curves', [0, 1, 2, 3], degfd=0)
    iess = ezt.define_essential(mesh, problem, 'curves', [0, 1, 2, 3], degfd=1,
                                iessp=iess)
    iess = ezt.define_essential(mesh, problem, 'points', [0], physq=1,
                                iessp=iess)

    # fill values for the essential boundary conditions

    uess = ezt.fill_system_vector(mesh, problem, 'curves', [2], func, funcnr=5)

    # apply essential boundary conditions to the system

    ezt.apply_essential(system_matrix, rhs, uess, iess)

    # solve the system

    u = spsolve(system_matrix.tocsr(), rhs)

    # pressure in all nodes for plotting

    xr = ezt.refcoor_nodal_points(mesh)
    user.psi, _ = ezt.basis_function('quad', 'Q1', xr)
    user.u = u
    derivs = {}
    derivs["pressure"] = ezt.deriv_vector(mesh, problem, ezt.stokes_pressure,
                                          user)

    # derivatives of the velocity

    xr = ezt.refcoor_nodal_points(mesh)
    user.phi, user.dphi = ezt.basis_function('quad', 'Q2', xr)
    user.u = u

    user.comp = 0  # dudx
    derivs["dudx"] = ezt.deriv_vector(mesh, problem, ezt.stokes_deriv, user)
    user.comp = 1  # dudy
    derivs["dudy"] = ezt.deriv_vector(mesh, problem, ezt.stokes_deriv, user)
    user.comp = 2  # dvdx
    derivs["dvdx"] = ezt.deriv_vector(mesh, problem, ezt.stokes_deriv, user)
    user.comp = 3  # dvdy
    derivs["dvdy"] = ezt.deriv_vector(mesh, problem, ezt.stokes_deriv, user)
    user.comp = 4  # dvdx - dudy = vorticity
    derivs["omega"] = ezt.deriv_vector(mesh, problem, ezt.stokes_deriv, user)
    user.comp = 6  # divu, divergence of the velocity field
    derivs["div_u"] = ezt.deriv_vector(mesh, problem, ezt.stokes_deriv, user)
    user.comp = 7  # gammadot, effective strain rate = sqrt(2II_D)
    derivs["gammadot"] = ezt.deriv_vector(mesh, problem, ezt.stokes_deriv,
                                          user)

    return max(derivs["omega"].u), mesh, problem, u, derivs


if __name__ == '__main__':
    results = main()
    print('max omega = ', results[0])
