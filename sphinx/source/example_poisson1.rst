example poisson1
================

.. code-block:: bash

    # Poisson problem on a unit square with Dirichlet boundary conditions

    import numpy as np
    import eztfem as ezt
    from func import func
    from scipy.sparse.linalg import spsolve

    # create mesh

    mesh = ezt.quadrilateral2d([20, 20], 'quad9')

    # define the problem

    elementdof = np.array([[1, 1, 1, 1, 1, 1, 1, 1, 1],
                        [2, 2, 2, 2, 2, 2, 2, 2, 2]], dtype=int).transpose()
    problem = ezt.Problem(mesh, elementdof, nphysq=1)

    # define Gauss integration and basis functions

    user = ezt.User()
    shape = 'quad'
    user.xr, user.wg = ezt.gauss_legendre(shape, n=3)
    user.phi, user.dphi = ezt.basis_function(shape, 'Q2', user.xr)

    # user struct for setting problem coefficients, ...

    user.coorsys = 0
    user.alpha = 1
    user.funcnr = 4
    user.func = func

    # assemble the system matrix and vector

    A, f = ezt.build_system(mesh, problem, ezt.poisson_elem, user)

    # define essential boundary conditions (Dirichlet)

    iess = ezt.define_essential(mesh, problem, 'curves', [0, 1, 2, 3])

    # fill values for the essential boundary conditions

    uess = ezt.fill_system_vector(mesh, problem, 'curves', [0, 1, 2, 3], func,
                                funcnr=3)

    # apply essential boundary conditions to the system

    ezt.apply_essential(A, f, uess, iess)

    # solve the system 

    u = spsolve(A.tocsr(), f)

    # compare with exact solution

    uex = ezt.fill_system_vector(mesh, problem, 'nodes',
                                np.arange(mesh.nnodes), func, funcnr=3)

    maxdiff = max(abs(u-uex))

    # gradient (dudx,dudy) of the solution 

    xr = ezt.refcoor_nodal_points(mesh)
    [user.phi, user.dphi] = ezt.basis_function('quad', 'Q2', xr)
    user.u = u
    gradu = ezt.deriv_vector(mesh, problem, ezt.poisson_deriv, user)