# Poisson problem on a square [-1,1]x[-1,1] with Dirichlet boundary conditions
# This problem solves Example 1.1.1 of the book of Elman et al.

import sys
sys.path.append('.')
sys.path.append('..')
sys.path.append('../..')

import numpy as np
import eztfem as ezt
from examples.poisson.func import func # use full path to be able to run from different folder
from scipy.sparse.linalg import spsolve

def main():

    # create mesh

    print('mesh')
    mesh = ezt.quadrilateral2d([40,40],'quad4',origin=np.array([-1,-1]),length=np.array([2,2]))

    # define the problem

    print('problem_definition')
    elementdof = np.array([[1,1,1,1],
                        [2,2,2,2]],dtype=int).transpose()
    problem = ezt.Problem(mesh,elementdof,nphysq=1)

    # define Gauss integration and basis functions

    user = ezt.User()
    shape='quad'

    print('gauss_legendre')
    user.xr, user.wg = ezt.gauss_legendre(shape,n=2)

    print('basis_function phi')
    user.phi, user.dphi = ezt.basis_function(shape,'Q1', user.xr )

    # user struct for setting problem coefficients, ...

    user.coorsys = 0
    user.alpha = 1
    user.funcnr = 5
    user.func = func

    # assemble the system matrix and vector

    print('build_system')
    A, f = ezt.build_system ( mesh, problem, ezt.poisson_elem, user )

    # define essential boundary conditions (Dirichlet)

    print('define_essential')
    iess = ezt.define_essential ( mesh, problem,'curves', [0,1,2,3] )

    # fill values for the essential boundary conditions

    print('fill_system_vector')
    uess = np.zeros([problem.numdegfd])

    # apply essential boundary conditions to the system

    print('apply_essential')
    ezt.apply_essential ( A, f, uess, iess )

    # solve the system 

    print('solve')
    u = spsolve(A.tocsr(), f)

    # maximum value

    print('Maximum value ' ,max(u))

    return max(u)

if __name__ == '__main__':
    main()