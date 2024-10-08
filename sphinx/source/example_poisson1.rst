Poisson problem with Dirichlet BCs
==================================

Importing packages
------------------
Import the relevant packages to run this script.

.. code-block:: bash

    # Poisson problem on a unit square with Dirichlet boundary conditions

    import numpy as np
    import eztfem as ezt
    from func import func
    from scipy.sparse.linalg import spsolve

Mesh generation
---------------
The mesh is created using a 20 by 20 mesh of nine-node quadrilaterals. 
The variable ``mesh`` is a structure having several 
components, including the coordinates of the nodes, the topology, the points 
and curves 
(see Figure XXXX). For more info on the ``Mesh``-structure 
see the package documentation and the documentation of the function 
``quadrilateral2d``. 

.. code-block:: bash

    # create mesh

    mesh = ezt.quadrilateral2d([20, 20], 'quad9')

Problem definition
------------------
The problem is defined. The variable ``problem`` is a structure having several 
components. The argument ``elementdof`` of the function ``problem_definition`` 
is a matrix where each column defines degrees of freedom in the nodes of an 
element. Note, that there is a transpose operator ``.transpose()`` in the 
assignment of ``elementdof``. The `physical quantities`` define the unknown
degrees of freedom in the system vector. In this case the number of physical
quantities is 1, which means that only the first column of ``elementdof`` 
defines the system vector of unknowns. The second column of ``elementdof`` is
only used for postprocessing (see below).

.. code-block:: bash

    # define the problem

    elementdof = np.array([[1, 1, 1, 1, 1, 1, 1, 1, 1],
                        [2, 2, 2, 2, 2, 2, 2, 2, 2]], dtype=int).transpose()
    problem = ezt.Problem(mesh, elementdof, nphysq=1)

Fill user object
----------------
The user object is used to pass various data and coefficients to the element
routines. Here, the Gauss-Legendre integration points (``xr``) and weights 
(``wg``) are computed and added, together with the basis functions 
:math:`\phi_k`, :math:`k=1,\dots,9`
and the derivatives :math:`\partial \phi_k/\partial \xi_j`,
:math:`k=1,\dots,9`, :math:`j=1,2` in these points for a bi-quadratic
(:math:`Q_2`) interpolation. The weights, basis functions and derivatives
are stored in the structure ``User`` for later use on element level.

.. code-block:: bash

    # define Gauss integration and basis functions

    user = ezt.User()
    shape = 'quad'
    user.xr, user.wg = ezt.gauss_legendre(shape, n=3)
    user.phi, user.dphi = ezt.basis_function(shape, 'Q2', user.xr)

Further information which is needed on elementlevel is added to ``user``:

* ``coorsys`` coordinate system: 0=Cartesian, 1=axisymmetric.
* ``alpha`` diffusion coefficient :math:`\alpha` which is 1 for the Poisson equation.
* ``funcnr`` function number in the function ``func`` for setting the right-hand side.
* ``func`` function handle to the function ``func``.

.. code-block:: bash

    # user struct for setting problem coefficients, ...

    user.coorsys = 0
    user.alpha = 1
    user.funcnr = 4
    user.func = func

System assembly
---------------
Assemble the system matrix :math:`\boldsymbol{A}` and vector 
:math:`\boldsymbol{f}` using the element function ``poisson_elem``.

.. code-block:: bash

    # assemble the system matrix and vector

    A, f = ezt.build_system(mesh, problem, ezt.poisson_elem, user)

Boundary conditions
-------------------
Define and apply Dirichlet boundary conditions. First, at line 32,
an index array ``iess`` is  
computed to indicate that the degrees ``u(iess)`` need to prescribed. Then,
the prescribed values 
are filled in the (system) vector ``uess``. Finally, the system matrix 
:math:`\boldsymbol{A}` and vector :math:`\boldsymbol{f}`
are modified to take the Dirichlet conditions into account.

.. code-block:: bash

    # define essential boundary conditions (Dirichlet)

    iess = ezt.define_essential(mesh, problem, 'curves', [0, 1, 2, 3])

    # fill values for the essential boundary conditions

    uess = ezt.fill_system_vector(mesh, problem, 'curves', [0, 1, 2, 3], func,
                                funcnr=3)

    # apply essential boundary conditions to the system

    ezt.apply_essential(A, f, uess, iess)

Solve linear system
-------------------
Solve the system :math:`\boldsymbol{A}\boldsymbol{u}=\boldsymbol{f}`.

.. code-block:: bash

    # solve the system 

    u = spsolve(A.tocsr(), f)

Compare solution
----------------
Print the maximum difference in the nodes, i.e. :math:`\max|u_i-u_{i,\text{exact}}|`.

.. code-block:: bash

    # compare with exact solution

    uex = ezt.fill_system_vector(mesh, problem, 'nodes',
                                np.arange(mesh.nnodes), func, funcnr=3)

    maxdiff = max(abs(u-uex))

Postprocessing
--------------
Derive a column vector (array) with :math:`\nabla u` in the nodes by averaging 
the values in elements
connected to the nodes. This column vector (array), ``gradu``, is defined by
the second column of ``elementdof`` and is a structure with two components: 
``gradu.vec`` the vector number 
(=column number in ``elementdof``) and ``gradu.u`` the actual data in all nodes. 
In order to derive
:math:`\vek\nabla u`, the basis functions and the derivatives of the basis 
functions need to be replaced by
the values in the nodes (using the function ``refcoor_nodal_points``). Also,
the system vector ``u`` needs to be available at the element level and is 
supplied via a component of ``user``.

.. code-block:: bash

    # gradient (dudx,dudy) of the solution 

    xr = ezt.refcoor_nodal_points(mesh)
    [user.phi, user.dphi] = ezt.basis_function('quad', 'Q2', xr)
    user.u = u
    gradu = ezt.deriv_vector(mesh, problem, ezt.poisson_deriv, user)