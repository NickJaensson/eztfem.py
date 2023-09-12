# Streamfunction from the velocity field

import sys
sys.path.append('/Users/njaensson/Desktop/eztfem.py/')

import numpy as np
from pprint import pprint
from src.quadrilateral2d import quadrilateral2d
from src.problem_class import Problem
from src.gauss_legendre import gauss_legendre
from src.basis_function import basis_function
from src.user_class import User
from examples.stokes.func import func

from src.build_system import build_system
from src.define_essential import define_essential
from src.fill_system_vector import fill_system_vector
from src.apply_essential import apply_essential
from src.refcoor_nodal_points import refcoor_nodal_points
from src.deriv_vector import deriv_vector

from addons.stokes.stokes_elem import stokes_elem
from addons.stokes.stokes_pressure import stokes_pressure
from addons.stokes.stokes_deriv import stokes_deriv
from addons.stokes.streamfunction_elem import streamfunction_elem

from scipy.sparse.linalg import spsolve

import pretty_errors

from src.pos_array import pos_array
from src.mesh_class import Mesh

from src.add_boundary_elements import add_boundary_elements
from addons.stokes.streamfunction_natboun_curve import streamfunction_natboun_curve

def main(mesh=None, problem=None, u=None):

    if mesh is None:
        raise ValueError("Error streamfunction1: mesh should be present")

    if problem is None:
        raise ValueError("Error streamfunction1: problem should be present")
    
    if u is None:
        raise ValueError("Error streamfunction1: u should be present")
    
    # define the problem

    print('problem_definition')
    elementdof = np.array([[1,1,1,1,1,1,1,1,1],
                           [2,2,2,2,2,2,2,2,2]],dtype=int).transpose()
    problem_s = Problem(mesh,elementdof,nphysq=1)

    # define Gauss integration and basis functions

    user = User()
    shape='quad'
    user.xr, user.wg = gauss_legendre(shape,n=3)
    user.phi, user.dphi = basis_function(shape,'Q2', user.xr )

    # user struct for setting problem coefficients, ...

    user.coorsys = 0

    # Get velocities node for node
    pos, _ = pos_array(problem, np.arange(mesh.nnodes), physq=0, order='ND')
    user.v = u[pos[0]]

    # assemble the system matrix and vector

    print('build_system')
    A, f = build_system ( mesh, problem_s, streamfunction_elem, user, posvectors=True )

    # define Gauss integration and basis functions (for boundary integral)

    print('gauss_legendre')
    [xr,user.wg]=gauss_legendre('line',n=3 )

    print('basis_function phi')
    [user.phi,user.dphi]=basis_function('line','P2', xr )

    # add natural boundary condition 

    print('add_boundary_elements')
    for crv in range(3):
        f = add_boundary_elements ( mesh, problem_s, f, 
                                   streamfunction_natboun_curve, user, 
                                   posvectors=True, curve=crv )
        
    # define essential boundary conditions (Dirichlet)

    print('define_essential')
    iess = define_essential ( mesh, problem_s, 'points',[0] )

    # fill values for the essential boundary conditions

    print('fill_system_vector')
    uess = np.zeros([problem_s.numdegfd])

    # apply essential boundary conditions to the system

    print('apply_essential')
    A, f, _ = apply_essential ( A, f, uess, iess )

    # solve the system 

    print('solve')
    streamf = spsolve(A.tocsr(), f)

    # maximum value (Matlab: 0.008522786557203 on 20x20 nelem)

    print('Maximum value ' ,max(streamf))

    return max(streamf) 

if __name__ == '__main__':

    import examples.stokes.stokes1

    _, mesh, problem, u = examples.stokes.stokes1.main()

    main(mesh,problem,u)