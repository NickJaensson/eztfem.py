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

from src_test.build_system import build_system
from src_test.poisson_elem import poisson_elem
from src_test.stokes_elem import stokes_elem

import pretty_errors

problemtype = "stokes"
mesh=quadrilateral2d([2,2],'quad9')

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
else:
    raise ValueError(f"Invalid problemtype : {problemtype}")

A = A.tocsr() # for similarity with TFEM
print(A)
print(f)