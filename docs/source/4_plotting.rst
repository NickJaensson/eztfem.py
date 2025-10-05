Chapter 3. Plotting
-------------------

As seen in the examples above, eztfem.py ships a set of basic plotting
functions. They open standard Matplotlib figures, so options like ``hold``
(``plt.figure``/``plt.show``), axes toggles, and styling can be customised
through normal Matplotlib APIs. The available plotting helpers include:

``plot_gauss_legendre``
    Plot Gauss–Legendre points on reference triangles and quadrilaterals.

``plot_basis_function``
    Plot basis functions on reference triangles and quadrilaterals.

``plot_points_curves``
    Plot all points and curves of the mesh with labelling included. This uses a
    special case of the more general ``plot_curves``.

``plot_curves``
    Plot curves, optionally restricted to a subset via keyword arguments.

``plot_mesh``
    Plot the mesh. Optional arguments allow plotting node markers, node
    numbers, and element numbers.

``plot_sol``
    Create a 2D or 3D (default) colour plot of the solution vector. Without any
    optional arguments the first degree of the first physical quantity is
    plotted.

``plot_sol_contour``
    Create a contour plot of the solution vector. The helper internally
    interpolates nodal data to a regular grid using ``scipy.interpolate``.

``plot_sol_quiver``
    Plot vector fields using quiver arrows. Supply the appropriate physical
    quantity (e.g. velocity components) to visualise flows.

``plot_sol_over_line``
    Plot a solution along a straight line. The sampled data can optionally be
    returned as NumPy arrays for further processing.

``plot_vector``
    Create a colour plot of a vector that was produced by ``deriv_vector``.
