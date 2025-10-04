Stokes problem with Neumann BCs
================================

We consider the Poisson problem on a rectangle. The problem is similar to 
``poisson1.m`` from the previous section, but now we change the domain to 
\((x, y) \in (0, 1.2) \times (0, 1)\)
and change the boundary condition on curve ``C2`` to a natural boundary. 
The latter means 
that we have to specify :math:`h_N = - \frac{\partial u}{\partial n}` on 
curve ``C2``.
We substitute the following exact solution for :math:`u`:

.. math::

   u = \cos(\pi x)\cos(\pi y) + x^3y^3

and find the corresponding :math:`f`:

.. math::

   f = 2\pi^2\cos(\pi x)\cos(\pi y) - 6(xy^3 + x^3y)

For the natural boundary condition we need:

.. math::

   h_\text{N} = - \frac{\partial u}{\partial n} = - \frac{\partial u}{\partial x} = \pi\sin(\pi x)\cos(\pi y) - 3x^2y^3

to be evaluated on curve ``C2``.

The Python script that solves this problem is ``poisson4.py`` in the 
``poisson`` folder of the examples folder. The main differences with respect
to \texttt{poisson1.m} will be discussed.

First, the essential boundary conditions need to be applied on 
curves 1, 3 and 4. Second, the boundary integral for the natural boundary 
condition needs be added:

.. code-block:: bash

    # define Gauss integration and basis functions (for boundary integral)

    print('gauss_legendre')
    [xr, user.wg] = ezt.gauss_legendre('line', n=3)

    print('basis_function phi')
    [user.phi, user.dphi] = ezt.basis_function('line', 'P2', xr)

    # add natural boundary condition

    print('add_boundary_elements')
    user.funcnr = 8
    ezt.add_boundary_elements(mesh, problem, f, ezt.poisson_natboun_curve,
                              user, curve=1)

Note, that the Gauss points and the basis functions and it's derivatives 
need to be redefined before using the boundary
element routine! Note also, that the right-hand side vector 
``f`` is part of the argument list for the function
``add_boundary_elements``.