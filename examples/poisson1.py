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

import pretty_errors

# create mesh

print('mesh')
mesh=quadrilateral2d([2,2],'quad9')

# define the problem

print('problem_definition')
elementdof=np.array([[1,1,1,1,1,1,1,1,1],
                     [2,2,2,2,2,2,2,2,2]],dtype=int).transpose()

problem=Problem(mesh,elementdof,nphysq=1)

# define Gauss integration and basis functions

print('gauss_legendre')
user = User()
user.xr, user.wg = gauss_legendre(shape='quad',n=3)
user.phi, user.dphi = basis_function('quad','Q2', user.xr )

# user object for setting problem coefficients, ...

user.coorsys = 0
user.alpha = 1
user.funcnr = 4
user.func = func

# assemble the system matrix and vector

print('build_system')
A, f = build_system ( mesh, problem, poisson_elem, user )
A = A.tocsr() # for similarity with TFEM

print(A)
print(f)