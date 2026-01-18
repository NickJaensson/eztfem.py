4. Classes
--------------------------

4.1 ``Mesh``
~~~~~~~~~~~~

The ``Mesh`` class contains:

``elnumnod``
    Number of nodes per element.
``elshape``
    Element shape identifier.
``eltype``
    Element type identifier.
``coords``
    Nodal coordinates.
``elnodes``
    Element-to-node connectivity array.
``nnodes``
    Number of nodes.
``nelem``
    Number of elements.
``points``
    Named points; useful for boundary conditions and data extraction.
``curves``
    Curve definitions that reference mesh nodes and line elements, enabling
    boundary conditions and data extraction.

Meshes can be generated with ``line1d`` (1D) and ``quadrilateral2d`` (2D).
Separate meshes can be merged into a single mesh with ``mesh_merge`` to create
more complicated geometries. Additional helpers are provided in the
``eztfem.addons.meshes`` module.

4.2 ``Problem``
~~~~~~~~~~~~~~~

The ``Problem`` structure stores data related to the degrees of freedom defined
on the mesh. Consult the documentation of ``eztfem.Problem`` for a complete
list. Some key components include:

``nvec``
    Number of vectors defined on the mesh. A vector is a storage container for
    degrees of freedom defined in the nodes and stored node-for-node.
``vec_elnumdegfd``
    Array giving the number of degrees of freedom in each node of an element
    for every vector. In Python this information matches the ``elementdof``
    array passed to ``Problem``.
``vec_nodnumdegfd``
    Number of degrees of freedom per vector in each node of the global mesh.
``vec_numdegfd``
    Number of degrees of freedom for each vector.
``nphysq``
    Number of physical degrees of freedom (i.e. those that are solved for in
    the system vector). The first ``nphysq`` columns of ``elementdof`` define
    the system vector; additional columns are available for post-processing.
``elnumdegfd``
    Number of system degrees of freedom in each nodal point of an element.
``numdegfd``
    Total number of system degrees of freedom.

4.3 ``Vector``
~~~~~~~~~~~
``eztfem.py`` distinguishes between *vectors* and *system vectors*. The system
vector contains all degrees of freedom that are solved for in the linear and
are of type NumPy array. Vectors are used for data storage and are of type
``Vector``. The ``Vector`` class is used for data storage, often as the result 
of ``deriv_vector``. A vector is
defined by one of the columns in ``elementdof`` (see Section 4.2) and . Vectors 
are storedas a small structure with two components (``vec`` and ``u``), where 
``vec`` is an index identifying the vector and ``u`` contains the actual data.

System vectors are defined by ``elnumdegfd`` and make up the degrees of freedom
that appear in the assembled linear system. They are stored as one-dimensional
NumPy arrays.

Use ``pos_array`` and ``pos_array_vec`` helpers to extract data from vectors and
system vectors. For example (assuming ``import numpy as np``)::

   idx = ezt.pos_array(problem, np.arange(3), physq=0, order='ND')

returns indices for the requested entries.
