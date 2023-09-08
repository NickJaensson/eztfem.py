# Poisson problem on an L-shape domain with Dirichlet boundary conditions
# This problem solves Example 1.1.2 of the book of Elman et al.

import sys
sys.path.append('/Users/njaensson/Desktop/eztfem.py/')

import numpy as np
from pprint import pprint
from src.quadrilateral2d import quadrilateral2d
from src.problem_class import Problem
from src.gauss_legendre import gauss_legendre
from src.basis_function import basis_function
from src.user_class import User
from examples.poisson.func import func

from src.build_system import build_system
from src.define_essential import define_essential
from src.fill_system_vector import fill_system_vector
from src.apply_essential import apply_essential
from src.refcoor_nodal_points import refcoor_nodal_points
from src.deriv_vector import deriv_vector

from addons.poisson.poisson_elem import poisson_elem
from addons.poisson.poisson_deriv import poisson_deriv

from scipy.sparse.linalg import spsolve

from addons.meshes.l_shape2d import l_shape2d

import pretty_errors

def main():

    # create mesh

    print('mesh')
    mesh = l_shape2d([20,20,20,20],'quad4',factor=np.array([10,10,10,10]))

    # define the problem

    print('problem_definition')
    elementdof = np.array([[1,1,1,1],
                        [2,2,2,2]],dtype=int).transpose()
    problem = Problem(mesh,elementdof,nphysq=1)

    # define Gauss integration and basis functions

    user = User()
    shape='quad'

    print('gauss_legendre')
    user.xr, user.wg = gauss_legendre(shape,n=2)

    print('basis_function phi')
    user.phi, user.dphi = basis_function(shape,'Q1', user.xr )

    # user struct for setting problem coefficients, ...

    user.coorsys = 0
    user.alpha = 1
    user.funcnr = 5
    user.func = func

    # assemble the system matrix and vector

    print('build_system')
    A, f = build_system ( mesh, problem, poisson_elem, user )

    # define essential boundary conditions (Dirichlet)

    print('define_essential')
    iess = define_essential ( mesh, problem,'curves', [0,1,2,3,4,5,6,7] )

    # fill values for the essential boundary conditions

    print('fill_system_vector')
    uess = np.zeros([problem.numdegfd])

    # apply essential boundary conditions to the system

    print('apply_essential')
    A, f, _ = apply_essential ( A, f, uess, iess )

    # solve the system 

    print('solve')
    u = spsolve(A.tocsr(), f)

    # maximum value

    print('Maximum value ' ,max(u))

    return max(u)

if __name__ == '__main__':
    main()