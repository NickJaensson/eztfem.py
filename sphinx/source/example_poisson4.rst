example poisson4
================

We consider the Poisson problem on a rectangle. The problem is similar to ``poisson1.m`` from the previous section, but now we change the domain to \((x, y) \in (0, 1.2) \times (0, 1)\)
and change the boundary condition on curve ``C2`` to a natural boundary. The latter means 
that we have to specify :math:`h_N = - \frac{\partial u}{\partial n}` on curve ``C2``.
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
