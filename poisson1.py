# Poisson problem on a unit square with Dirichlet boundary conditions
# Manufactured solution.

#clear all; close all;

from quadrilateral2d import quadrilateral2d
from distribute_elements import distribute_elements
from problem_definition import Problem
import numpy as np

from pprint import pprint

# create mesh

print('mesh')
mesh=quadrilateral2d([1,2],'quad9')

# define the problem

print('problem_definition')
elementdof=np.array([[1,1,1,1,1,1,1,1,1],
                     [2,2,2,2,2,2,2,2,2]],dtype=int)

problem=Problem(mesh,elementdof,nphysq=1)

pprint(vars(problem))



