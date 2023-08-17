# Poisson problem on a unit square with Dirichlet boundary conditions
# Manufactured solution.

import sys
sys.path.append('../')

import numpy as np
from pprint import pprint
from src.quadrilateral2d import quadrilateral2d
from src.problem_definition import Problem
from src.gauss_legendre import gauss_legendre
from src.basis_function import basis_function
from src.user import User
from func import func

from src.build_system import build_system
from src.define_essential import define_essential
from src.fill_system_vector import fill_system_vector
from src.apply_essential import apply_essential
from src.refcoor_nodal_points import refcoor_nodal_points
from src.deriv_vector import deriv_vector

from addons.poisson.poisson_elem import poisson_elem
from addons.stokes.stokes_elem import stokes_elem
from addons.stokes.stokes_pressure import stokes_pressure
from addons.stokes.stokes_deriv import stokes_deriv

from scipy.sparse.linalg import spsolve

import pretty_errors

problemtype = "stokes"
mesh=quadrilateral2d([3,2],'quad9')

if problemtype == "poisson":
    elementdof=np.array([[1,1,1,1,1,1,1,1,1],
                         [2,2,2,2,2,2,2,2,2]],dtype=int).transpose()
    problem=Problem(mesh,elementdof,nphysq=1)

    user = User()
    shape='quad'
    user.xr, user.wg = gauss_legendre(shape,n=3)
    user.phi, user.dphi = basis_function(shape,'Q2', user.xr )

    user.coorsys = 0
    user.alpha = 1
    user.funcnr = 4
    user.func = func

    A, f = build_system ( mesh, problem, poisson_elem, user )

    iess = define_essential ( mesh, problem,'curves', [0,1,2,3] )

    uess = fill_system_vector ( mesh, problem, 'curves', [0,1], func, funcnr=3 )
    uess = fill_system_vector ( mesh, problem, 'curves', [2,3], func, funcnr=3, fin=uess )

    A, f, _ = apply_essential ( A, f, uess, iess )

    u = spsolve(A.tocsr(), f)


elif problemtype == "stokes":
    elementdof=np.array([[2,2,2,2,2,2,2,2,2],
                         [1,0,1,0,1,0,1,0,0],
                         [1,1,1,1,1,1,1,1,1]],dtype=int).transpose()
    problem=Problem(mesh,elementdof,nphysq=2)

    user = User()
    shape='quad'
    user.xr, user.wg = gauss_legendre(shape,n=3)
    user.phi, user.dphi = basis_function(shape,'Q2', user.xr )
    user.psi, _ = basis_function(shape,'Q1', user.xr )

    user.coorsys = 0
    user.mu = 1
    user.funcnr = 0 
    user.func = func # not used when funcnr == 0

    A, f = build_system ( mesh, problem, stokes_elem, user )

    iess = define_essential ( mesh, problem, 'curves',[0,1,2,3], degfd=0 )
    iess = define_essential ( mesh, problem, 'curves',[0,2,2,3], degfd=1, iessp=iess )
    iess = define_essential ( mesh, problem, 'points',[0], physq=1, iessp=iess )

    uess = fill_system_vector ( mesh, problem, 'curves', [0,1], func, funcnr=3 )
    uess = fill_system_vector ( mesh, problem, 'curves', [2,3], func, funcnr=3, fin=uess )

    A, f, _ = apply_essential ( A, f, uess, iess )

    u = spsolve(A.tocsr(), f)

    xr = refcoor_nodal_points ( mesh )
    user.psi, _ = basis_function('quad','Q1', xr )
    user.u = u
    pressure = deriv_vector ( mesh, problem, stokes_pressure, user )

    user.phi, user.dphi = basis_function('quad','Q2', xr )
    user.comp = 6 # divu, divergence of the velocity field
    divu = deriv_vector ( mesh, problem, stokes_deriv, user ) 
    user.comp = 7 # gammadot, effective strain rate = sqrt(2II_D) 
    gammadot = deriv_vector ( mesh, problem, stokes_deriv, user )

else:
    raise ValueError(f"Invalid problemtype : {problemtype}")