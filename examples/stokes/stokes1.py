# Poisson problem on a unit square with Dirichlet boundary conditions
# Manufactured solution.

import sys
sys.path.append('/Users/njaensson/Desktop/eztfem.py/')

import numpy as np
from pprint import pprint
from src.quadrilateral2d import quadrilateral2d
from src.problem_class import Problem
from src.gauss_legendre import gauss_legendre
from src.basis_function import basis_function
from src.user_class import User
from examples.stokes.func import func

from src.build_system import build_system
from src.define_essential import define_essential
from src.fill_system_vector import fill_system_vector
from src.apply_essential import apply_essential
from src.refcoor_nodal_points import refcoor_nodal_points
from src.deriv_vector import deriv_vector

from addons.stokes.stokes_elem import stokes_elem
from addons.stokes.stokes_pressure import stokes_pressure
from addons.stokes.stokes_deriv import stokes_deriv

from scipy.sparse.linalg import spsolve

import pretty_errors

def main():

    # create mesh

    print('mesh')
    mesh = quadrilateral2d([20,20],'quad9')

    # define the problem

    print('problem_definition')
    elementdof = np.array([[2,2,2,2,2,2,2,2,2],
                        [1,0,1,0,1,0,1,0,0],
                        [1,1,1,1,1,1,1,1,1]],dtype=int).transpose()
    problem = Problem(mesh,elementdof,nphysq=2)

    # define Gauss integration and basis functions

    user = User()
    shape='quad'
    user.xr, user.wg = gauss_legendre(shape,n=3)
    user.phi, user.dphi = basis_function(shape,'Q2', user.xr )
    user.psi, _ = basis_function(shape,'Q1', user.xr )

    # user struct for setting problem coefficients, ...

    user.coorsys = 0
    user.mu = 1
    user.funcnr = 0 
    user.func = func # not used when funcnr == 0

    # assemble the system matrix and vector

    print('build_system')
    A, f = build_system ( mesh, problem, stokes_elem, user )

    # define essential boundary conditions (Dirichlet)

    print('define_essential')
    iess = define_essential ( mesh, problem, 'curves',[0,1,2,3], degfd=0 )
    iess = define_essential ( mesh, problem, 'curves',[0,2,2,3], degfd=1, iessp=iess )
    iess = define_essential ( mesh, problem, 'points',[0], physq=1, iessp=iess )

    # fill values for the essential boundary conditions

    print('fill_system_vector')
    uess = fill_system_vector ( mesh, problem, 'curves', [0,1], func, funcnr=3 )
    uess = fill_system_vector ( mesh, problem, 'curves', [2,3], func, funcnr=3, fin=uess )

    # apply essential boundary conditions to the system

    print('apply_essential')
    A, f, _ = apply_essential ( A, f, uess, iess )

    # solve the system 

    print('solve')
    u = spsolve(A.tocsr(), f)

    # Pressure in all nodes for plotting

    print('pressure in nodes')
    xr = refcoor_nodal_points ( mesh )
    user.psi, _ = basis_function('quad','Q1', xr )
    user.u = u
    pressure = deriv_vector ( mesh, problem, stokes_pressure, user )

    # derivatives of the velocity

    print('velocity derivatives')
    xr = refcoor_nodal_points ( mesh )
    user.phi, user.dphi = basis_function('quad','Q2', xr )
    user.u = u

    user.comp = 0 # dudx
    dudx = deriv_vector ( mesh, problem, stokes_deriv, user ) 
    user.comp = 1 # dudy
    dudy = deriv_vector ( mesh, problem, stokes_deriv, user )
    user.comp = 2 # dvdx
    dvdx = deriv_vector ( mesh, problem, stokes_deriv, user ) 
    user.comp = 3 # dvdy
    dvdy = deriv_vector ( mesh, problem, stokes_deriv, user )
    user.comp = 4 # dvdx - dudy = vorticity
    omega = deriv_vector ( mesh, problem, stokes_deriv, user ) 
    user.comp = 6 # divu, divergence of the velocity field
    divu = deriv_vector ( mesh, problem, stokes_deriv, user ) 
    user.comp = 7 # gammadot, effective strain rate = sqrt(2II_D) 
    gammadot = deriv_vector ( mesh, problem, stokes_deriv, user )

    return max(gammadot.u)

if __name__ == '__main__':
    main()