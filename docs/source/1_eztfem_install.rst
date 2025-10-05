eztfem installation
===================

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


In the file browser you see the following folders in the eztfem tree:

* ``src``
* ``examples``

The folder ``src`` contains the folders ``core`` and ``addons``.
The folder ``core`` contains core functions of eztfem, whereas the folder
``addons`` contains the add-ons, with a separate folder for each add-on. 
The following add-ons are available: 

* ``plotlib``
* ``poisson``
* ``stokes``
* ``meshes``

which contain plot functions, element routines for the Poisson and Stokes 
equations and functions for more complicated meshes, respectively. 
The folder ``examples`` contain the examples, divided into separate 
folders. Each folder might contain several examples, but all of the same 
kind, e.g. a problem involving the Stokes equation. The following example 
folders are available: 

* ``poisson``
* ``stokes``
