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

from scipy.sparse.linalg import spsolve

import pretty_errors

mesh = quadrilateral2d([3,2],'quad9')

elementdof = np.array([[1,1,1,1,1,1,1,1,1],
                       [2,2,2,2,2,2,2,2,2]],dtype=int).transpose()
problem = Problem(mesh,elementdof,nphysq=1)

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