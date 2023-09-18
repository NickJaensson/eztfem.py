# Poisson problem on a unit square with Dirichlet boundary conditions
# Manufactured solution.

import sys
sys.path.append('..')

import numpy as np
import eztfem as ezt
from examples.stokes.func import func # use full path to be able to run from different folder
from scipy.sparse.linalg import spsolve

def main():

    # create mesh

    print('mesh')
    mesh = ezt.quadrilateral2d([20,20],'quad9')

    # define the problem

    print('problem_definition')
    elementdof = np.array([[2,2,2,2,2,2,2,2,2],
                        [1,0,1,0,1,0,1,0,0],
                        [1,1,1,1,1,1,1,1,1]],dtype=int).transpose()
    problem = ezt.Problem(mesh,elementdof,nphysq=2)

    # define Gauss integration and basis functions

    user = ezt.User()
    shape='quad'

    print('gauss_legendre')
    user.xr, user.wg = ezt.gauss_legendre(shape,n=3)

    print('basis_function phi')
    user.phi, user.dphi = ezt.basis_function(shape,'Q2', user.xr )

    print('basis_function psi')
    user.psi, _ = ezt.basis_function(shape,'Q1', user.xr )

    # user struct for setting problem coefficients, ...

    user.coorsys = 0
    user.mu = 1
    user.funcnr = 0 
    user.func = func # not used when funcnr == 0

    # assemble the system matrix and vector

    print('build_system')
    A, f = ezt.build_system ( mesh, problem, ezt.stokes_elem, user )

    # define essential boundary conditions (Dirichlet)

    print('define_essential')
    iess = ezt.define_essential ( mesh, problem, 'curves',[0,1,2,3], degfd=0 )
    iess = ezt.define_essential ( mesh, problem, 'curves',[0,1,2,3], degfd=1, iessp=iess )
    iess = ezt.define_essential ( mesh, problem, 'points',[0], physq=1, iessp=iess )

    # fill values for the essential boundary conditions

    print('fill_system_vector')
    uess = ezt.fill_system_vector ( mesh, problem, 'curves', [2], func, funcnr=5 )

    # apply essential boundary conditions to the system

    print('apply_essential')
    A, f, _ = ezt.apply_essential ( A, f, uess, iess )

    # solve the system 

    print('solve')
    u = spsolve(A.tocsr(), f)

    # Pressure in all nodes for plotting

    print('pressure in nodes')
    xr = ezt.refcoor_nodal_points ( mesh )
    user.psi, _ = ezt.basis_function('quad','Q1', xr )
    user.u = u
    pressure = ezt.deriv_vector ( mesh, problem, ezt.stokes_pressure, user )

    # derivatives of the velocity

    print('velocity derivatives')
    xr = ezt.refcoor_nodal_points ( mesh )
    user.phi, user.dphi = ezt.basis_function('quad','Q2', xr )
    user.u = u

    user.comp = 0 # dudx
    dudx = ezt.deriv_vector ( mesh, problem, ezt.stokes_deriv, user ) 
    user.comp = 1 # dudy
    dudy = ezt.deriv_vector ( mesh, problem, ezt.stokes_deriv, user )
    user.comp = 2 # dvdx
    dvdx = ezt.deriv_vector ( mesh, problem, ezt.stokes_deriv, user ) 
    user.comp = 3 # dvdy
    dvdy = ezt.deriv_vector ( mesh, problem, ezt.stokes_deriv, user )
    user.comp = 4 # dvdx - dudy = vorticity
    omega = ezt.deriv_vector ( mesh, problem, ezt.stokes_deriv, user ) 
    user.comp = 6 # divu, divergence of the velocity field
    divu = ezt.deriv_vector ( mesh, problem, ezt.stokes_deriv, user ) 
    user.comp = 7 # gammadot, effective strain rate = sqrt(2II_D) 
    gammadot = ezt.deriv_vector ( mesh, problem, ezt.stokes_deriv, user )

    return max(omega.u), mesh, problem, u

if __name__ == '__main__':
    main()