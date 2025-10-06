1. Introduction
-----------------------

1.1 What is eztfem.py?
~~~~~~~~~~~~~~~~~~~~~~

Eztfem.py is a simple toolkit for the finite element method (FEM) intended for
use in teaching finite elements for fluid flow. It can be seen as a
significantly reduced version of TFEM [1]. Since eztfem.py consists of Python
modules only, it is easily accessible for students.

1.2 Installation
~~~~~~~~~~~~~~~~

Install the eztfem package by running the following command in the eztfem.py
root directory:

.. code-block:: bash

   pip install -e .

Note: this will install the pacakge in `editable` mode, which means that the
package links directly to the code (instead of copying it to the Python 
environment), thus changes to the package will be reflected immediately.

All functions in eztfem can be found in these docs. There are two ways of 
calling eztfem functions. 

* The first option imports the entire package and calls functions by

.. code-block:: bash

   import eztfem as ezt
   ezt.FUNCTION_NAME()

* The second option import specific functions and call them by

.. code-block:: bash

   from eztfem import FUNCTION_NAME
   FUNCTION_NAME()

Note that in the second option, all functions need to be imported explicitly.

1.3 Folder layout of eztfem.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

At the top level of the repository you will find, among others, the following
folders:

``eztfem``
    The package folder containing the core implementation. Within this folder
    the ``core`` subpackage contains the core functionality of eztfem.py,
    whereas the ``addons`` subpackage contains add-ons, each in its own module.
    The following add-ons are available: ``plotlib`` (plotting helpers),
    ``poisson`` and ``stokes`` (problem-specific element routines), and
    ``meshes`` (utilities for more complicated meshes).

``examples``
    Contains the example problems divided into separate folders. Each folder
    can contain several examples, but all of the same kind, e.g. problems
    involving the Stokes equation. The following example folders are available:
    ``poisson`` and ``stokes``.

``docs``
    Sphinx documentation sources and configuration.

1.4 Running examples
~~~~~~~~~~~~~~~~~~~~

All examples are available as normal Python scripts (``.py``) and as 
Jupyter Notebooks (``.ipynb``). For post-processing (e.g., plotting), we 
recommend using the Jupyter Notebooks, since data will remain available 
after running functions.

Using VS Code, Jupyter Notebooks can be run directly after installing the 
relevant extensions. Python scripts can be run by clicking the Play button 
above the script. This will open a terminal and run the script with the Python 
version as defined by your VS Code setup. 

For more details on the examples, see Section ``examples``.

1.5 Documentation
~~~~~~~~~~~~~~~~~

Documentation on the functions in eztfem.py can be obtained using the built-in
Python ``help`` function or via ``pydoc``. For example, from the Python REPL:

.. code-block:: pycon

   >>> import eztfem as ezt
   >>> help(ezt.build_system)
