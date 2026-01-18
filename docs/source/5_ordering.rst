5. Vector and matrix ordering
--------------------------

5.1 Working with degrees of freedom
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In `eztfem` the full system vector and vectors are stored in *NPD* ordering, 
i.e., nodes are in the outer loop, physical quantities in the 
middle loop and degrees of freedom in the inner loop. For example, if the
velocity is the first physical quantity (with degrees of freedom :math:`u` and 
:math:`v` and pressure is the second physical quantity (with degree of freedom 
:math:`p`), the system vector would be stored as

.. math::

   [u_1, v_1, p_1, u_2, v_2, p_2, u_3, v_3, p_3]

When building the system matrix, it has some advantages to use a different 
ordering on element level. In `eztfem`, the default in `build_system` assumes
that elements are written in *PDN* ordering, which is given by

.. math::

    [u_1, u_2, u_3, v_1, v_2, v_3, p_1, p_2, p_3]

Note, that this only applies to how the element routines are written, i.e.,
it does not change the ordering of the full system vector itself 
(which is always in *NPD* ordering).

5.2 System matrix and right-hand side
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

eztfem.py imposes Dirichlet conditions following the approach of [2]. Let the
full system of equations be partitioned into unknown (:math:`u`) and prescribed
(:math:`p`) degrees of freedom:

.. math::

   \begin{bmatrix} A_{uu} & A_{up} \\ A_{pu} & A_{pp} \end{bmatrix}
   \begin{bmatrix} \tilde{u}_u \\ \tilde{u}_p \end{bmatrix}
   =
   \begin{bmatrix} \tilde{f}_u \\ \tilde{f}_p \end{bmatrix}.

After substituting the Dirichlet data :math:`\tilde{u}_p = \tilde{u}_D` we
obtain the reduced system

.. math::

   A_{uu} \tilde{u}_u = \tilde{f}_u - A_{up} \tilde{u}_D.

Instead of assembling and solving the reduced system directly, eztfem.py
modifies the full system as

.. math::

   \begin{bmatrix} A_{uu} & 0 \\ 0 & I \end{bmatrix}
   \begin{bmatrix} \tilde{u}_u \\ \tilde{u}_p \end{bmatrix}
   =
   \begin{bmatrix} \tilde{f}_u - A_{up} \tilde{u}_D \\ \tilde{u}_D \end{bmatrix},

where :math:`I` is the identity matrix corresponding to the prescribed degrees
of freedom. This is implemented in ``apply_essential``. Although the resulting
system is slightly larger, it is straightforward to implement. Be aware that
this approach adds eigenvalues equal to 1 with multiplicity equal to the number
of prescribed degrees of freedom. Note, that in practice, the system matrix is
not reordered in ``apply_essential``, i.e., the matrix is modified in place,
keeping the original *NPD* ordering.
