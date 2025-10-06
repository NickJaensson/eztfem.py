2. Examples
-------------------

This chapter explains several example problems in more detail. A typical
eztfem.py script for a problem is structured as follows:

.. figure:: figs/script.pdf
   :alt: alt text
   :align: center
   :width: 60%

   **Figure 2.1:** Typical structure of an eztfem.py script for solving a 
   problem.


2.1 A Poisson problem on a unit square with Dirichlet boundary conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Consider a Poisson problem:

.. math::

   -\nabla^2 u = f

on a unit square :math:`(x,y) \in (0,1) \times (0,1)`. Substituting the
solution

.. math::

   u = 1 + \cos(\pi x)\cos(\pi y)

gives the corresponding right-hand side:

.. math::

   f = 2\pi^2 \cos(\pi x)\cos(\pi y)

The script that solves this problem is ``poisson1.py`` in the ``poisson``
folder of the ``examples`` directory. The script is listed below and explained
line by line after the listing.

.. literalinclude:: ../../examples/poisson/poisson1.py
   :language: python
   :linenos:

Lines 15
    The mesh is created using a 20×20 grid of nine-node quadrilaterals. The
    variable ``mesh`` is a structure with several components, including the
    coordinates of the nodes, the topology, the points and the curves. For more
    information on the mesh structure see Section 4.1 and the documentation of
    ``quadrilateral2d``. The function ``quadrilateral2d`` has several options
    to create a mesh on a quadrilateral.

.. figure:: figs/unitsquare.pdf
   :alt: alt text
   :align: center
   :width: 60%

   **Figure 2.2:** Points and curves defined by the mesh generator function
   ``quadrilateral2d``. Points and curves are basically only placeholders for
   certain nodes and element edges respectively. Note the default directions of
   the curves.

Lines 20–22
    The problem is defined. The variable ``problem`` is a structure with
    several components. The argument ``elementdof`` is a matrix where each
    column defines degrees of freedom in the nodes of an element. Note that the
    call to ``np.array(...).transpose()`` makes the column-major layout explicit
    so that the first column defines the system vector of unknowns. The second
    column of ``elementdof`` is only used for post-processing (see below).

Lines 26–33
    A ``User`` instance is filled with Gauss-Legendre integration points
    (``user.xr``) and weights (``user.wg``), together with the basis functions
    :math:`\phi_k`, :math:`k=1,\ldots,9` and the derivatives
    :math:`\partial \phi_k/\partial \xi_j`, :math:`k=1,\ldots,9`,
    :math:`j=1,2` in these points for a bi-quadratic (:math:`Q_2`)
    interpolation. The weights, basis functions and derivatives are stored in
    the ``user`` structure for later use on element level.

Lines 37–40
    The ``user`` structure is filled with additional data needed at element
    level: ``coorsys`` (coordinate system: 0 = Cartesian, 1 = axisymmetric),
    ``alpha`` (diffusion coefficient :math:`\alpha`, which is 1 for the Poisson
    equation), ``funcnr`` (function selector for the right-hand side), and
    ``func`` (callable providing the forcing terms).

Lines 45
    The system matrix :math:`\boldsymbol{A}` and vector :math:`\boldsymbol{f}`
    are assembled using ``eztfem.build_system`` with the Poisson element
    routine.

Lines 50-55
    Dirichlet boundary conditions are defined on all boundary curves and the
    prescribed values are filled into ``uess`` using ``fill_system_vector``.
    ``apply_essential`` modifies :math:`\boldsymbol{A}` and :math:`\boldsymbol{f}`
    to impose the essential boundary conditions.

Lines 61
    The linear system is solved with SciPy’s ``spsolve``. The maximum
    difference between the numerical and analytical solution in the nodes is
    printed.

Lines 80-84
    The gradient :math:`(\partial u/\partial x, \partial u/\partial y)` is
    computed. The basis functions are recomputed in the nodal points and the
    solution ``u`` is supplied to the element routine via the ``user``
    structure. The gradient is returned in ``gradu``.

To visualise the results you can use the plotting helpers, for example:

- Plot the mesh:

.. code-block:: python

    ezt.plot_mesh(mesh)

- Plot the solution :math:`u`:

.. code-block:: python

    ezt.plot_sol(mesh, problem, u)

- Plot the solution :math:`u` over a line (e.g., from point [0,0] to point 
  [1,1]):

.. code-block:: python

    ezt.plot_sol_over_line(mesh, problem, u, [[0.0, 0.0], [1.0, 1.0]])

- Plot the first component of the gradient (:math:`\partial u / \partial x`):

.. code-block:: python

    ezt.plot_vector(mesh, problem, gradu, degfd=0)


2.2 A Poisson problem on a rectangle with Dirichlet and Neumann boundary conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This problem is similar to ``poisson1.py`` but with a natural boundary 
condition on one side. The exact solution is

.. math::

   u = \cos(\pi x)\cos(\pi y) + x^3 y^3

which yields the forcing term

.. math::

   f = 2\pi^2 \cos(\pi x)\cos(\pi y) - 6(xy^3 + x^3y)

The natural boundary condition on curve ``C1`` (see Figure 2.2) requires

.. math::

   h_N = -\frac{\partial u}{\partial n} = -\frac{\partial u}{\partial x}
       = \pi \sin(\pi x)\cos(\pi y) - 3x^2 y^3

which must be evaluated along that curve.

The corresponding script is ``poisson4.py``. Compared with ``poisson1.py`` the
main differences are:

* Essential boundary conditions are imposed on curves 0, 2 and 3 only (curve
  indices are zero-based in eztfem.py).
* ``user.funcnr`` values 7 and 8 are used for the manufactured forcing term and
  the Neumann boundary contribution, respectively.
* Before adding the boundary integral for the natural boundary condition, the
  Gauss points and basis functions for the boundary element must be defined:

  .. code-block:: python

     user.xr, user.wg = ezt.gauss_legendre('line', n=3)
     user.phi, user.dphi = ezt.basis_function('line', 'P2', user.xr)

* The natural boundary contribution is added with
  ``ezt.add_boundary_elements`` on curve 1 (the second curve in zero-based
  indexing). Note that ``f`` is part of the argument list::

     user.funcnr = 8
     ezt.add_boundary_elements(
         mesh, problem, f, ezt.poisson_natboun_curve, user, curve=1
     )

2.3 A Stokes problem on a unit square with Dirichlet boundary conditions: driven cavity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We consider the steady Stokes equations

.. math::

   -\nabla \cdot (2\mu\boldsymbol{D}) + \nabla p = \boldsymbol{f}, \qquad
   \nabla \cdot \boldsymbol{u} = 0,

where :math:`\boldsymbol{D} = (\nabla \boldsymbol{u} + \nabla \boldsymbol{u}^T)/2`,
for a driven cavity problem on the unit square. 

.. figure:: figs/lid_driven.pdf
   :alt: alt text
   :align: center
   :width: 60%

   **Figure 2.3:** Driven cavity problem.

On the upper boundary curve ``C2`` (see Figure 2.2) requires a
horizontal velocity of 1 is imposed and on all other walls the velocity
components are zero. In the lower-left corner the pressure value is set to zero
(Dirichlet).

The script ``stokes1.py`` in ``examples/stokes`` solves this problem. Important
changes with respect to the Poisson example include:

* The problem definition uses two physical quantities (velocity and pressure):

  .. code-block:: python

     elementdof = np.array([
         [2, 2, 2, 2, 2, 2, 2, 2, 2],
         [1, 0, 1, 0, 1, 0, 1, 0, 0],
         [1, 1, 1, 1, 1, 1, 1, 1, 1],
     ], dtype=int).transpose()
     problem = ezt.Problem(mesh, elementdof, nphysq=2)

* Velocity boundary conditions are prescribed with ``fill_system_vector`` on
  curve 2 (the top boundary) and combined with the pressure constraint at point
  0 (handled through ``define_essential``) before calling ``apply_essential``.
* Plotting requires supplying the ``problem`` definition so that vector fields
  can be visualised, e.g. ``ezt.plot_sol_quiver(mesh, problem, u)``.

2.4 Flow generated by a traction (pressure) difference: Poiseuille
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This problem is similar to the driven-cavity setup but introduces traction
boundary conditions (see ``stokes2.py``) to generate a Poiseuille flow.

.. figure:: figs/poiseuille.pdf
   :alt: alt text
   :align: center
   :width: 60%

   **Figure 2.4:** Poiseuille flow problem.

The main differences compared with
``stokes1.py`` are:

* A traction boundary condition is added on curve``C3`` using
  ``ezt.add_boundary_elements`` with the Stokes boundary element routine and
  ``physqrow=[0]``.
* ``user.func`` is switched to ``traction_func`` (with ``funcnr = 1``) to
  supply the traction data.
* Dirichlet conditions for the velocity components are imposed using separate
  calls to ``define_essential``: :math:`u` on curves 0 and 2, and
  :math:`v` on curves ``C0`` to ``C3``.
* The flow rate through curve 1 (the second curve) is computed with
  ``ezt.integrate_boundary_elements`` and the element routine
  ``stokes_flowrate_curve``. The basis functions on the curve are the same as
  those used when adding the traction contribution.

2.5 Computing the streamfunction field from a known velocity field
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For plotting streamlines in 2D flows we compute the streamfunction
:math:`\psi` from the velocity vector 
:math:`\boldsymbol{u} = u\boldsymbol{e}_1 + v\boldsymbol{e}_2`.
The definition of :math:`\psi` is

.. math::

   \frac{\partial \psi}{\partial x} = -v, \qquad
   \frac{\partial \psi}{\partial y} = u,

which shows that

.. math::

   \psi(x_2) - \psi(x_1)
   = \int_{\Gamma_{12}} \boldsymbol{n} \cdot \boldsymbol{u} \, ds
   \qquad

or in words: the streamfunction difference is the flow rate through a curve
connecting two points. The equation above can be used to compute :math:`\psi` 
in all nodal points, but the procedure is cumbersome and not straightforward.
Therefore, we prefer a procedure that solves the Poisson equation

.. math::

   -\nabla^2 \psi = \omega \qquad

with the vorticity 
:math:`\omega = \partial v/\partial x - \partial u/\partial y`.
We use a Neumann boundary condition:

.. math::

   \frac{\partial \psi}{\partial n}
   = \boldsymbol{n} \cdot \nabla \psi
   = \boldsymbol{n} \cdot (-v\boldsymbol{e}_1 + u\boldsymbol{e}_2)
   = -v n_1 + u n_2 = -\boldsymbol{u} \cdot \boldsymbol{t} \qquad

where :math:`\boldsymbol{t}` is the tangential vector 
(with :math:`\boldsymbol{n}` to the ‘right’).

For plotting streamlines in an axisymmetrical flow we have

.. math::

   \frac{\partial \psi}{\partial z} = -2\pi r v, \qquad
   \frac{\partial \psi}{\partial r} = 2\pi r u \qquad

where :math:`(z, r)` are cylindrical coordinates and :math:`u` is the velocity
component in the axial (:math:`z`) direction and :math:`v` the velocity
component in the radial (:math:`r`) direction
(:math:`\boldsymbol{u} = u\boldsymbol{e}_z + v\boldsymbol{e}_r`).
We easily find that

.. math::

   \psi(x_2) - \psi(x_1)
   = \int_{\Gamma_{12}} \boldsymbol{n} \cdot \boldsymbol{v} \, 2\pi r \, ds

or in words: the streamfunction difference is the flow rate through a surface
consisting of a curve connecting two points that is expanded over the full
circumferential (:math:`\theta`) direction. We find from the equations that

.. math::

   -\nabla^2_{zr} \psi
   = -\left(
       \frac{\partial^2 \psi}{\partial z^2}
       + \frac{\partial^2 \psi}{\partial r^2}
     \right)
   = 2\pi r \left(\frac{\partial v}{\partial z} 
     - \frac{\partial u}{\partial r}\right)
   = 2\pi (r\omega + u)

with the vorticity :math:`\omega = \partial v/\partial z - \partial u/\partial r`.
Note that :math:`\nabla^2_{zr}` is the *planar* (:math:`z, r`) Laplace operator
and not the real :math:`\nabla^2` expressed in a cylindrical coordinate system.
We use a Neumann boundary condition:

.. math::

   \frac{\partial \psi}{\partial n}
   = \boldsymbol{n} \cdot \nabla \psi
   = 2\pi r \boldsymbol{n} \cdot (-v\boldsymbol{e}_z + u\boldsymbol{e}_r)
   = 2\pi r(-v n_z + u n_r)
   = -2\pi r \boldsymbol{u} \cdot \boldsymbol{t}
   \qquad

where :math:`\boldsymbol{t}` is the tangential vector 
(with :math:`\boldsymbol{n}` to the ‘right’).
Note that the Neumann condition is zero on the center line (:math:`r = 0`).



See ``streamfunction1.py`` for the implementation details and remember that
``stokes1.py`` must be executed first to provide the velocity field.

Important eztfem.py features highlighted by ``streamfunction1.py``:

* A new ``Problem`` instance is created for the scalar streamfunction while a
  second column in ``elementdof`` is used to import the velocity field from the
  Stokes solution.
* ``pos_array`` obtains the positions of the velocity solution in the system
  vector so that the values can be supplied through the ``user`` structure.
* ``build_system`` receives the optional keyword argument ``posvectors`` to
  make vector data available at element level.
* ``add_boundary_elements`` assembles Neumann contributions on each curve.
