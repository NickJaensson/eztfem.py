# Poisson problem on a unit square with Dirichlet boundary conditions
# Manufactured solution.

from quadrilateral2d import quadrilateral2d
from distribute_elements import distribute_elements
from problem_definition import Problem
from gauss_legendre import gauss_legendre
from basis_function import basis_function
import numpy as np
from user import User

from pprint import pprint

# create mesh

print('mesh')
mesh=quadrilateral2d([1,2],'quad9')
#pprint(vars(mesh))

# define the problem

print('problem_definition')
elementdof=np.array([[1,1,1,1,1,1,1,1,1],
                     [2,2,2,2,2,2,2,2,2]],dtype=int).transpose()

problem=Problem(mesh,elementdof,nphysq=1)

#pprint(vars(problem))

# define Gauss integration and basis functions

print('gauss_legendre')
user = User()
user.xr, user.wg = gauss_legendre(shape='quad',n=3)

# Usage example
user.phi, user.dphi = basis_function('quad','Q2', user.xr )

print(user.phi,user.dphi)