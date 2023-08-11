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

# create mesh

print('mesh')
mesh=quadrilateral2d([1,2],'quad9')

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