# Streamfunction from the velocity field

import sys
sys.path.append('..')

import numpy as np
import eztfem as ezt
from scipy.sparse.linalg import spsolve

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
    problem_s = ezt.Problem(mesh,elementdof,nphysq=1)

    # define Gauss integration and basis functions

    user = ezt.User()
    shape='quad'
    user.xr, user.wg = ezt.gauss_legendre(shape,n=3)
    user.phi, user.dphi = ezt.basis_function(shape,'Q2', user.xr )

    # user struct for setting problem coefficients, ...

    user.coorsys = 0

    # Get velocities node for node
    pos, _ = ezt.pos_array(problem, np.arange(mesh.nnodes), physq=0, order='ND')
    user.v = u[pos[0]]

    # assemble the system matrix and vector

    print('build_system')
    A, f = ezt.build_system ( mesh, problem_s, ezt.streamfunction_elem, user, posvectors=True )

    # define Gauss integration and basis functions (for boundary integral)

    print('gauss_legendre')
    [xr,user.wg]=ezt.gauss_legendre('line',n=3 )

    print('basis_function phi')
    [user.phi,user.dphi]=ezt.basis_function('line','P2', xr )

    # add natural boundary condition 

    print('add_boundary_elements')
    for crv in range(3):
        ezt.add_boundary_elements ( mesh, problem_s, f, 
                                   ezt.streamfunction_natboun_curve, user, 
                                   posvectors=True, curve=crv )
        
    # define essential boundary conditions (Dirichlet)

    print('define_essential')
    iess = ezt.define_essential ( mesh, problem_s, 'points',[0] )

    # fill values for the essential boundary conditions

    print('fill_system_vector')
    uess = np.zeros([problem_s.numdegfd])

    # apply essential boundary conditions to the system

    print('apply_essential')
    ezt.apply_essential ( A, f, uess, iess )

    # solve the system 

    print('solve')
    streamf = spsolve(A.tocsr(), f)

    # maximum value (Matlab: 0.008522786557203 on 20x20 nelem)

    print('Maximum value ' ,max(streamf))

    return max(streamf) 

if __name__ == '__main__':

    import stokes1

    _, mesh, problem, u = stokes1.main()

    main(mesh,problem,u)