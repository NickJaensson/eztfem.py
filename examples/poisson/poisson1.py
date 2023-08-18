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
from addons.poisson.poisson_deriv import poisson_deriv

from scipy.sparse.linalg import spsolve

import pretty_errors

# create mesh

print('mesh')
mesh = quadrilateral2d([2,2],'quad9')

# define the problem

print('problem_definition')
elementdof = np.array([[1,1,1,1,1,1,1,1,1],
                       [2,2,2,2,2,2,2,2,2]],dtype=int).transpose()
problem = Problem(mesh,elementdof,nphysq=1)

# define Gauss integration and basis functions

user = User()
shape='quad'

print('gauss_legendre')
user.xr, user.wg = gauss_legendre(shape,n=3)

print('basis_function phi')
user.phi, user.dphi = basis_function(shape,'Q2', user.xr )

# user struct for setting problem coefficients, ...

user.coorsys = 0
user.alpha = 1
user.funcnr = 4
user.func = func

# assemble the system matrix and vector

print('build_system')
A, f = build_system ( mesh, problem, poisson_elem, user )

# define essential boundary conditions (Dirichlet)

print('define_essential')
iess = define_essential ( mesh, problem,'curves', [0,1,2,3] )

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

# compare with exact solution

print('Difference with exact solution:')
uex = fill_system_vector ( mesh, problem, 'nodes', np.arange(mesh.nnodes), func, funcnr=3 ) 

print(max(abs(u-uex)))

# gradient (dudx,dudy) of the solution 

print('gradient (dudx,dudy)')
xr = refcoor_nodal_points ( mesh )
[user.phi,user.dphi]=basis_function('quad','Q2', xr )
user.u = u

gradu = deriv_vector ( mesh, problem, poisson_deriv, user )
print(gradu.u)