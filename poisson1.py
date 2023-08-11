# Poisson problem on a unit square with Dirichlet boundary conditions
# Manufactured solution.

from quadrilateral2d import quadrilateral2d
from distribute_elements import distribute_elements
from problem_definition import Problem
from gauss_legendre import gauss_legendre
import numpy as np

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
xr, wg = gauss_legendre(shape='quad',n=3)

#print(xr)
#print(wg)

