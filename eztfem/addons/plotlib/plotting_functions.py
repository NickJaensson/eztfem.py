'''
Module for plotting functions using PyVista.
'''
import copy
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
from pyvista import CellType
from ...core.pos_array import pos_array, pos_array_vec
from ...core.shapefunc import basis_function
from ...core.gauss import gauss_legendre


def fill_mesh_pv(mesh_pv, problem, u, physq, degfd):
    """
    Fills the point data of a mesh object with the values from a given solution
    array based on the specified degrees of freedom.

    Parameters
    ----------
    mesh_pv : pyvista.PolyData
        A PyVista mesh

    problem : Problem
        Eztfem problem object

    u : array_like
        A NumPy array representing the solution values for the degrees of
        freedom (DOFs) that will be mapped to the mesh.

    physq : int
        Physical quantity identifier, which should be an integer.

    degfd : int, list, or np.ndarray
        Degree of freedom indices to be plotted. If a single integer is
        provided, the function maps the values for that DOF to the mesh. If a
        list or array of two integers is provided, it maps two DOFs to the
        mesh.

    Raises
    ------
    ValueError
        If `degfd` is not of length 1 or 2.

    Notes
    -----
    A deep copy of the object is made to avoid modifying the mesh_pv argument.

    """

    mesh_pv_plot = copy.deepcopy(mesh_pv)

    assert isinstance(physq, (int, np.integer))

    if isinstance(degfd, (int, np.integer)):
        degfd = [degfd]

    if isinstance(degfd, np.ndarray):
        degfd = degfd.tolist()

    nnodes = mesh_pv_plot.number_of_points

    match len(degfd):
        case 1:
            u_plot = np.zeros(nnodes)

            for node in range(nnodes):
                posn, _ = pos_array(problem, node, physq=physq, order='DN')
                u_plot[node] = u[posn[0][degfd[0]]]

            mesh_pv_plot.point_data['u'] = u_plot

        case 2:
            u0 = np.zeros(nnodes)
            u1 = np.zeros(nnodes)

            for node in range(nnodes):
                posn, _ = pos_array(problem, node, physq=physq, order='DN')
                u0[node] = u[posn[0][0]]
                u1[node] = u[posn[0][1]]

            mesh_pv_plot.point_data['u'] = \
                np.transpose(np.vstack((u0, u1, np.zeros_like(u0))))

        case _:
            raise ValueError("degfd must be of length 1 (int) or 2")

    return mesh_pv_plot


def fill_mesh_pv_vector(mesh_pv, problem, vector, degfd):
    """
    Fills the point data of a mesh object with the values from a given solution
    array based on the specified degrees of freedom.

    Parameters
    ----------
    mesh_pv : pyvista.PolyData
        A PyVista mesh

    problem : Problem
        Eztfem problem object

    vector : Vector
        The derived type vector to fill.

    degfd : int
        Degree of freedom indices to be added to mesh_pv

    Notes
    -----
    A deep copy of the object is made to avoid modifying the mesh_pv argument.

    """

    mesh_pv_plot = copy.deepcopy(mesh_pv)

    assert isinstance(degfd, (int, np.integer))

    nnodes = mesh_pv_plot.number_of_points

    u_plot = np.zeros(nnodes)

    for node in range(nnodes):
        posn, _ = pos_array_vec(problem, node, vec=vector.vec, order='DN')
        u_plot[node] = vector.u[posn[0][degfd]]

    mesh_pv_plot.point_data['u'] = u_plot

    return mesh_pv_plot


def plot_mesh_pv(mesh_pv, **kwargs):
    """
    Plot a mesh.

    Parameters
    ----------
    mesh_pv : pyvista.PolyData
        The mesh to be plotted.

    Keyword arguments
    -----------------
    kwargs : dict, optional
        Additional keyword arguments to pass to the plotter.add_mesh function.

    Notes
    -----
    This tutorial was followed to hide the internal edges:
    https://github.com/pyvista/pyvista/discussions/5777

    """

    style = kwargs.get('style', "wireframe")
    color = kwargs.get('color', "black")
    window_size = kwargs.get('window_size', (800, 400))

    kwargs.pop('style', None)
    kwargs.pop('color', None)
    kwargs.pop('window_size', None)

    surface = mesh_pv.separate_cells().extract_surface(nonlinear_subdivision=4)
    edges = surface.extract_feature_edges()

    plotter = pv.Plotter(window_size=window_size)
    plotter.add_mesh(surface)
    plotter.add_mesh(edges, style=style, color=color, **kwargs)
    plotter.camera_position = 'xy'
    plotter.show()


def plot_sol(mesh_pv, problem, u, **kwargs):
    """
    Plots the solution of a given problem on a mesh using PyVista.

    Parameters
    ----------
    mesh_pv : pyvista.PolyData
        The mesh on which to plot the solution.
    problem : Problem
        The problem object.
    u : numpy.ndarray
        The solution vector.

    Keyword arguments
    -----------------
    kwargs : dict, optional
        Additional keyword arguments to pass to the plotter.add_mesh function.

    """

    # Optional arguments
    physq = kwargs.get('physq', 0)
    degfd = kwargs.get('degfd', 0)
    window_size = kwargs.get('window_size', (800, 400))

    kwargs.pop('physq', None)
    kwargs.pop('degfd', None)
    kwargs.pop('window_size', None)

    mesh_pv_plot = fill_mesh_pv(mesh_pv, problem, u, physq, degfd)

    plotter = pv.Plotter(window_size=window_size)
    plotter.add_mesh(mesh_pv_plot, scalars="u", **kwargs)
    plotter.camera_position = 'xy'
    plotter.add_text((f'sol physq = {physq:d}  degfd = {physq:d}'),
                     font_size=12)
    plotter.show()


def plot_sol_contour(mesh_pv, problem, u, nlevels=10, **kwargs):
    """
    Plots the solution of a given problem on a mesh using PyVista.

    Parameters
    ----------
    mesh_pv : pyvista.PolyData
        The mesh on which to plot the solution.
    problem : Problem
        The problem object.
    u : numpy.ndarray
        The solution vector.
    nlevels : int
        Number of contour levels to plot (default = 10)

    Keyword arguments
    -----------------
    kwargs : dict, optional
        Additional keyword arguments to pass to the plotter.add_mesh function.

    """

    # Optional arguments
    physq = kwargs.get('physq', 0)
    degfd = kwargs.get('degfd', 0)
    window_size = kwargs.get('window_size', (800, 400))

    kwargs.pop('physq', None)
    kwargs.pop('degfd', None)
    kwargs.pop('window_size', None)

    mesh_pv_plot = fill_mesh_pv(mesh_pv, problem, u, physq, degfd)

    contours = mesh_pv_plot.contour(nlevels, scalars='u')

    plotter = pv.Plotter(window_size=window_size)
    plotter.add_mesh(mesh_pv_plot, color="lightgrey", **kwargs)
    plotter.add_mesh(contours)  # color="black", line_width=1)
    plotter.camera_position = 'xy'
    plotter.add_text((f'sol physq = {physq:d}  degfd = {physq:d}'),
                     font_size=12)
    plotter.show()


def plot_vector(mesh_pv, problem, vector, degfd=0, **kwargs):
    """
    Plots the solution of a given problem on a mesh using PyVista.

    Parameters
    ----------
    mesh_pv : pyvista.PolyData
        The mesh on which to plot the solution.
    problem : Problem
        The problem object.
    u : Vector
        The vector object.
    degfd : int
        The degree of freedom to plot (default = 0)

    Keyword arguments
    -----------------
    kwargs : dict, optional
        Additional keyword arguments to pass to the plotter.add_mesh function.

    """

    # Optional arguments
    degfd = kwargs.get('degfd', 0)
    window_size = kwargs.get('window_size', (800, 400))

    kwargs.pop('degfd', None)
    kwargs.pop('window_size', None)

    mesh_pv_plot = fill_mesh_pv_vector(mesh_pv, problem, vector, degfd)

    plotter = pv.Plotter(window_size=window_size)
    plotter.add_mesh(mesh_pv_plot, scalars="u", **kwargs)
    plotter.camera_position = 'xy'
    plotter.add_text((f'sol vec = {vector.vec:d}  degfd = {degfd:d}'),
                     font_size=12)
    plotter.show()


def plot_vector_contours(mesh_pv, problem, vector, degfd=0, nlevels=10,
                         **kwargs):
    """
    Plots the contours of a solution of a given problem on a mesh using
    PyVista.

    Parameters
    ----------
    mesh_pv : pyvista.PolyData
        The mesh on which to plot the solution.
    problem : Problem
        The problem object.
    u : Vector
        The vector object.
    degfd : int
        The degree of freedom to plot (default = 0)
    nlevels : int
        Number of contour levels to plot (default = 10)

    Keyword arguments
    -----------------
    kwargs : dict, optional
        Additional keyword arguments to pass to the plotter.add_mesh function.

    """

    # Optional arguments
    degfd = kwargs.get('degfd', 0)
    window_size = kwargs.get('window_size', (800, 400))

    kwargs.pop('degfd', None)
    kwargs.pop('window_size', None)

    mesh_pv_plot = fill_mesh_pv_vector(mesh_pv, problem, vector, degfd)

    contours = mesh_pv_plot.contour(nlevels, scalars='u')

    plotter = pv.Plotter(window_size=window_size)
    plotter.add_mesh(mesh_pv_plot, color="lightgrey", **kwargs)
    plotter.add_mesh(contours)  # color="black", line_width=1)
    plotter.camera_position = 'xy'
    plotter.add_text((f'sol vec = {vector.vec:d}  degfd = {degfd:d}'),
                     font_size=12)
    plotter.show()


def plot_mesh_plt(mesh, **kwargs):
    """
    Plot mesh structure.

    Parameters:
    mesh: mesh structure
    Optional arguments:
    - nodemarks: plot node marks (default: 0)
    - nodenumbers: plot node numbers (default: 0)
    - elementnumbers: plot element numbers (default: 0)
    """

    if mesh.ndim != 2:
        raise ValueError('Only 2D meshes can be plotted')

    # Set default options
    nodemarks = kwargs.get('nodemarks', 0)
    nodenumbers = kwargs.get('nodenumbers', 0)
    elementnumbers = kwargs.get('elementnumbers', 0)

    # Determine element type
    if mesh.elshape in [3, 4, 5]:  # all boundary nodes
        elnumnod = mesh.elnumnod
    elif mesh.elshape in [6, 7, 9, 10]:  # do not include center nodes
        elnumnod = mesh.elnumnod - 1
    else:
        raise ValueError(f'Invalid elshape = {mesh.elshape}')

    # Set coordinates
    x = np.zeros((elnumnod, mesh.nelem))
    y = np.zeros((elnumnod, mesh.nelem))

    for elem in range(mesh.nelem):
        for node in range(elnumnod):
            x[node, elem] = mesh.coor[mesh.topology[node, elem], 0]
            y[node, elem] = mesh.coor[mesh.topology[node, elem], 1]

    # Plot figure
    plt.fill(x, y, facecolor='none', edgecolor='k')
    plt.axis('equal')
    plt.axis('off')
    plt.title('Mesh')

    # Plot node marks
    if nodemarks:
        for node in range(mesh.nnodes):
            plt.text(mesh.coor[node, 0], mesh.coor[node, 1], 'X', fontsize=12,
                     color='g', ha='center')

    # Plot node numbers
    if nodenumbers:
        for node in range(mesh.nnodes):
            plt.text(mesh.coor[node, 0], mesh.coor[node, 1], str(node),
                     fontsize=8, color='b')

    # Plot element numbers
    if elementnumbers:
        for elem in range(mesh.nelem):
            plt.text(np.mean(x[:, elem]), np.mean(y[:, elem]), str(elem),
                     fontsize=8, color='r', ha='center', va='center')

    plt.show()


def plot_curves(mesh, **kwargs):
    """
    Plot curves (and optionally points) of the mesh.

    Parameters:
    mesh: mesh structure
    Optional arguments:
    - curves: the curves to be plotted (default: all curves)
    - nodemarks: plot node marks (default: 0)
    - nodenumbers: plot node numbers (0=off, 1=local, 2=global; default: 0)
    - elementnumbers: plot element numbers (default: 0)
    - curvenumbers: plot curve numbers (default: 0)
    - pointnumbers: plot point numbers (default: 0)
    """

    if mesh.ndim != 2:
        raise ValueError('Only 2D curves can be plotted')

    # Set default options
    curves = list(range(mesh.ncurves))
    nodemarks = 0
    nodenumbers = 0
    elementnumbers = 0
    curvenumbers = 0
    pointnumbers = 0

    # Update options based on kwargs
    curves = kwargs.get('curves', curves)
    nodemarks = kwargs.get('nodemarks', nodemarks)
    nodenumbers = kwargs.get('nodenumbers', nodenumbers)
    elementnumbers = kwargs.get('elementnumbers', elementnumbers)
    curvenumbers = kwargs.get('curvenumbers', curvenumbers)
    pointnumbers = kwargs.get('pointnumbers', pointnumbers)

    # Set coordinates (all curves must have the same element type)
    elnumnod = mesh.curves[0].elnumnod
    nelem = 0
    for curve in curves:
        if elnumnod != mesh.curves[curve].elnumnod:
            raise ValueError('All curves must have the same element type')
        nelem += mesh.curves[curve].nelem

    x = np.zeros((elnumnod, nelem))
    y = np.zeros((elnumnod, nelem))

    el = 0
    for curve in curves:
        for elem in range(mesh.curves[curve].nelem):
            el += 1
            for node in range(elnumnod):
                x[node, el-1] = \
                    mesh.coor[mesh.curves[curve].topology[node, elem, 1], 0]
                y[node, el-1] = \
                    mesh.coor[mesh.curves[curve].topology[node, elem, 1], 1]

    # Plot figure
    plt.plot(x, y, color='k')
    plt.axis('equal')
    plt.axis('off')

    # Plot node marks
    if nodemarks:
        for curve in curves:
            for node in mesh.curves[curve].nodes:
                plt.text(mesh.coor[node, 0], mesh.coor[node, 1], 'X',
                         fontsize=12, color='g', ha='center')

    # Plot node numbers
    if nodenumbers == 1:
        for curve in curves:
            for nod, node in enumerate(mesh.curves[curve].nodes):
                plt.text(mesh.coor[node, 0], mesh.coor[node, 1], str(nod + 1),
                         fontsize=8, color='b')
    elif nodenumbers == 2:
        for curve in curves:
            for node in mesh.curves[curve].nodes:
                plt.text(mesh.coor[node, 0], mesh.coor[node, 1], str(node + 1),
                         fontsize=8, color='b')

    # Plot element numbers
    if elementnumbers:
        el = 0
        for curve in curves:
            for elem in range(mesh.curves[curve].nelem):
                el += 1
                plt.text(np.mean(x[:, el-1]), np.mean(y[:, el-1]),
                         str(elem + 1), fontsize=8, color='r', ha='center',
                         va='center')

    # Plot curve numbers
    if curvenumbers:
        for curve in curves:
            x1 = np.min(mesh.coor[mesh.curves[curve].nodes, :], axis=0)
            x2 = np.max(mesh.coor[mesh.curves[curve].nodes, :], axis=0)
            plt.text((x1[0] + x2[0]) / 2, (x1[1] + x2[1]) / 2,
                     'C' + str(curve), fontsize=8, color='r',
                     ha='center', va='center')

    # Plot point numbers
    if pointnumbers:
        for point in range(mesh.npoints):
            node = mesh.points[point]
            plt.text(mesh.coor[node, 0], mesh.coor[node, 1],
                     'P' + str(point), fontsize=8, color='b')

    plt.show()


def plot_points_curves(mesh):
    """Plot all points and curves .

    Parameters
    ----------
    mesh : Mesh
        Mesh object containing the curves and points to be plotted.
    """

    plot_curves(mesh, pointnumbers=1, curvenumbers=1)


def plot_sol_over_line(mesh_pv, problem, u, points, physq=0, degfd=0,
                       npoints=200, plot_mesh=False):
    """
    Plots and samples data along a line through a mesh, optionally visualizing
    the mesh and the line.

    Parameters
    ----------
    mesh_pv : object
        A PyVista mesh object containing the data to be sampled and visualized.

    points : list or array_like
        A list of two 3D points, each of shape (3,), defining the start and end
        points of the line over which to sample.

    npoints : int, optional
        The number of points to sample along the line. Default is 200.

    plot_mesh : bool, optional
        If True, the function plots the mesh along with the line for
        visualization. The default is False.

    Returns
    -------
    sampled_data : object
        A PyVista dataset containing the sampled data along the line.

    """

    assert isinstance(physq, int)
    assert isinstance(degfd, int)

    if plot_mesh:
        line = pv.Line(points[0], points[1])
        p = pv.Plotter(window_size=(800, 400))
        p.add_mesh(mesh_pv, color="w")
        p.add_mesh(line, color="b")
        p.camera_position = 'xy'
        p.show()

    mesh_pv_plot = fill_mesh_pv(mesh_pv, problem, u, physq, degfd)

    mesh_pv_plot.plot_over_line(points[0], points[1], resolution=npoints)

    return mesh_pv_plot.sample_over_line(points[0], points[1],
                                         resolution=npoints)


def plot_vector_over_line(mesh_pv, problem, vector, points, degfd=0,
                          npoints=200, plot_mesh=False):
    """
    Plots and samples data along a line through a mesh, optionally visualizing
    the mesh and the line.

    Parameters
    ----------
    mesh_pv : object
        A PyVista mesh object containing the data to be sampled and visualized.

    points : list or array_like
        A list of two 3D points, each of shape (3,), defining the start and end
        points of the line over which to sample.

    npoints : int, optional
        The number of points to sample along the line. Default is 200.

    plot_mesh : bool, optional
        If True, the function plots the mesh along with the line for
        visualization. The default is False.

    Returns
    -------
    sampled_data : object
        A PyVista dataset containing the sampled data along the line.

    """

    assert isinstance(degfd, int)

    if plot_mesh:
        line = pv.Line(points[0], points[1])
        p = pv.Plotter(window_size=(800, 400))
        p.add_mesh(mesh_pv, color="w")
        p.add_mesh(line, color="b")
        p.camera_position = 'xy'
        p.show()

    mesh_pv_plot = fill_mesh_pv_vector(mesh_pv, problem, vector, degfd)

    mesh_pv_plot.plot_over_line(points[0], points[1], resolution=npoints)

    return mesh_pv_plot.sample_over_line(points[0], points[1],
                                         resolution=npoints)


def plot_quiver(mesh_pv, problem, u, **kwargs):
    """
    Plots the solution of a given problem on a mesh using PyVista.

    Parameters
    ----------
    mesh_pv : pyvista.PolyData
        The mesh on which to plot the solution.
    problem : Problem
        The problem object.
    u : numpy.ndarray
        The solution vector.
    scale : float optional
        Scaling factor for the arrows (default = 0.1)

    Keyword arguments
    -----------------
    kwargs : dict, optional
        Additional keyword arguments to pass to the plotter.add_mesh function.

    """

    # Optional arguments
    physq = kwargs.get('physq', 0)
    window_size = kwargs.get('window_size', (800, 400))
    scale = kwargs.get('scale', 0.1)

    kwargs.pop('physq', None)
    kwargs.pop('degfd', None)
    kwargs.pop('window_size', None)

    ndf = problem.elementdof[0, physq]

    assert (problem.elementdof[0, physq] == 2)
    assert np.all(problem.elementdof[0, physq] == ndf)

    mesh_pv_plot = fill_mesh_pv(mesh_pv, problem, u, physq, degfd=[0, 1])

    glyphs = mesh_pv_plot.glyph(orient="u", scale=True, factor=scale)

    plotter = pv.Plotter(window_size=window_size)
    # plotter.add_mesh(glyphs, show_scalar_bar=False, lighting=False,
    #                  cmap='coolwarm')
    plotter.add_mesh(mesh_pv_plot, color="lightgrey")
    plotter.add_mesh(glyphs, color="black")
    plotter.camera_position = 'xy'
    plotter.add_text((f'sol physq = {physq:d}  degfd = {physq:d}'),
                     font_size=12)
    plotter.show()


def plot_basis_function(shape, intpol, degfd, *, n=10, plot3d=True,
                        edges=True, axisoff=False, show=True, **kwargs):
    """Plot a reference-element basis function using PyVista.

    Parameters
    ----------
    shape : {'quad', 'triangle'}
        Shape of the reference element.
    intpol : {'P0', 'P1', 'P1+', 'P2', 'P2+', 'Q1', 'Q1+', 'Q2'}
        Interpolation family used to compute the basis functions.
    degfd : int
        Degree of freedom index of the basis function to plot. Both zero-based
        and one-based indexing are supported for convenience.

    Keyword Arguments
    -----------------
    n : int, optional
        Number of subdivisions in each coordinate direction. Default ``10``.
    plot3d : bool, optional
        Plot the basis function as a 3D surface (default) or as a coloured 2D
        patch when set to ``False``.
    edges : bool, optional
        If ``True`` (default) display the sub-element edges.
    axisoff : bool, optional
        Hide axes and scalar bar when ``True``. Default ``False``.
    show : bool, optional
        If ``True`` (default) immediately render the scene. Set to ``False`` to
        obtain the :class:`pyvista.Plotter` and mesh for further customization.
    window_size : tuple, optional
        Size of the PyVista rendering window. Defaults to ``(800, 400)``.
    kwargs : dict, optional
        Additional keyword arguments forwarded to
        :meth:`pyvista.Plotter.add_mesh`.

    Returns
    -------
    plotter : pyvista.Plotter
        Plotter containing the rendered basis function.
    mesh : pyvista.UnstructuredGrid
        Generated mesh with point data ``'phi'`` holding the basis values.

    """

    if not isinstance(n, (int, np.integer)) or n <= 0:
        raise ValueError("n must be a positive integer")

    if shape not in {"quad", "triangle"}:
        raise ValueError(f"Invalid shape: {shape}")

    if not isinstance(degfd, (int, np.integer)):
        raise TypeError("degfd must be an integer")

    if shape == "quad":
        points2d, cells, celltypes = _reference_quad_mesh(n)
    else:
        points2d, cells, celltypes = _reference_triangle_mesh(n)

    phi, _ = basis_function(shape, intpol, points2d)

    if degfd < 0:
        raise ValueError("degfd must be non-negative")

    if degfd >= phi.shape[1]:
        if degfd == phi.shape[1]:
            degfd_idx = degfd - 1
        else:
            raise ValueError(
                f"degfd={degfd} outside valid range [0, {phi.shape[1] - 1}]"
            )
    else:
        degfd_idx = int(degfd)

    phi_values = phi[:, degfd_idx]

    if plot3d:
        z = phi_values
    else:
        z = np.zeros_like(phi_values)

    points = np.column_stack((points2d, z))
    mesh = pv.UnstructuredGrid(cells, celltypes, points)
    mesh.point_data['phi'] = phi_values

    window_size = kwargs.pop('window_size', (800, 400))
    plotter = pv.Plotter(window_size=window_size)

    mesh_kwargs = {
        'scalars': 'phi',
        'show_edges': edges,
    }
    mesh_kwargs.update(kwargs)

    plotter.add_mesh(mesh, **mesh_kwargs)

    if not plot3d:
        plotter.view_xy()

    title = f"phi      node = {degfd_idx + 1}   intpol = {intpol}"
    plotter.add_title(title)

    if axisoff:
        plotter.remove_scalar_bar()
    else:
        plotter.show_bounds()

    if show:
        plotter.show()

    return plotter, mesh


def _reference_quad_mesh(n):
    """Create reference coordinates and topology for a subdivided quad."""

    nn1 = n + 1
    delta = 2.0 / n

    points = np.zeros((nn1 * nn1, 2))

    for j in range(nn1):
        for i in range(nn1):
            node = j * nn1 + i
            points[node, 0] = -1.0 + i * delta
            points[node, 1] = -1.0 + j * delta

    cells = []

    for j in range(n):
        row = j * nn1
        next_row = (j + 1) * nn1
        for i in range(n):
            n0 = row + i
            n1 = row + i + 1
            n2 = next_row + i + 1
            n3 = next_row + i
            cells.extend([4, n0, n1, n2, n3])

    cell_array = np.array(cells, dtype=np.int64)
    celltypes = np.full(n * n, CellType.QUAD, dtype=np.uint8)

    return points, cell_array, celltypes


def _reference_triangle_mesh(n):
    """Create reference coordinates and topology for a subdivided triangle."""

    delta = 1.0 / n
    nn1 = n + 1

    nnodes = nn1 * (nn1 + 1) // 2
    points = np.zeros((nnodes, 2))

    node = 0
    for j in range(nn1):
        for i in range(nn1 - j):
            points[node, 0] = i * delta
            points[node, 1] = j * delta
            node += 1

    cells = []
    celltypes = []

    np_idx = 0
    for j in range(n):
        first_count = n - j
        for i in range(first_count):
            n0 = np_idx + i
            n1 = np_idx + i + 1
            n2 = np_idx + i + n - j + 1
            cells.extend([3, n0, n1, n2])
            celltypes.append(CellType.TRIANGLE)

        second_count = n - j - 1
        for i in range(second_count):
            n0 = np_idx + i + 1
            n1 = np_idx + i + n - j + 2
            n2 = np_idx + i + n - j + 1
            cells.extend([3, n0, n1, n2])
            celltypes.append(CellType.TRIANGLE)

        np_idx += n - j + 1

    cell_array = np.array(cells, dtype=np.int64)
    celltypes = np.array(celltypes, dtype=np.uint8)

    return points, cell_array, celltypes


def plot_gauss_legendre(shape, **kwargs):
    """Plot Gauss-Legendre integration points.

    Parameters
    ----------
    shape : {'quad', 'triangle'}
        Shape of the domain. ``'quad'`` corresponds to the square
        :math:`[-1, 1] \times [-1, 1]` and ``'triangle'`` corresponds to the
        reference triangle defined by the lower-left half of
        :math:`[0, 1] \times [0, 1]`.

    Keyword Arguments
    -----------------
    n : int, optional
        Number of integration points per direction. Required when
        ``shape`` is ``'quad'``.
    p : int, optional
        Order of the integration rule. Required when ``shape`` is
        ``'triangle'``.
    ax : matplotlib.axes.Axes, optional
        Existing axes to plot on. When omitted, a new figure and axes are
        created.
    marker : str, optional
        Marker used to display the integration points (default: ``'+'``).
    color : str, optional
        Marker color (default: ``'k'``).
    markersize : float, optional
        Size of the marker (default: ``8``).
    show : bool, optional
        Display the plot immediately (default: ``True``).
    kwargs
        Additional keyword arguments forwarded
        to :func:`matplotlib.axes.Axes.plot`.

    Returns
    -------
    matplotlib.axes.Axes
        Axes containing the plot of the integration points.
    """

    n = kwargs.pop('n', None)
    p = kwargs.pop('p', None)
    ax = kwargs.pop('ax', None)
    marker = kwargs.pop('marker', '+')
    color = kwargs.pop('color', 'k')
    markersize = kwargs.pop('markersize', 8)
    show = kwargs.pop('show', True)

    if not isinstance(shape, str):
        raise ValueError('shape must be a string')

    shape = shape.lower()

    if shape == 'quad':
        if n is None:
            raise ValueError('n must be specified for the integration rule')
        points, _ = gauss_legendre('quad', n=n)
        x_limits = (-1, 1)
        y_limits = (-1, 1)
    elif shape == 'triangle':
        if p is None:
            raise ValueError('p must be specified for the integration rule')
        points, _ = gauss_legendre('triangle', p=p)
        x_limits = (0, 1)
        y_limits = (0, 1)
    else:
        raise ValueError(f'Invalid shape: {shape}')

    if ax is None:
        _, ax = plt.subplots()

    plot_kwargs = kwargs
    plot_kwargs.setdefault('linestyle', 'None')

    ax.plot(points[:, 0], points[:, 1], marker=marker, color=color,
            markersize=markersize, **plot_kwargs)

    ax.set_xlim(*x_limits)
    ax.set_ylim(*y_limits)

    if shape == 'triangle':
        ax.plot([1, 0], [0, 1], color='k')

    ax.set_aspect('equal', adjustable='box')

    if show:
        plt.show()

    return ax
