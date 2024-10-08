import numpy as np
from ...core.pos_array import pos_array
import pyvista as pv
import matplotlib.pyplot as plt
import copy


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


def plot_mesh(mesh, **kwargs):
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
            plt.text(mesh.coor[node, 0], mesh.coor[node, 1], str(node + 1),
                     fontsize=8, color='b')

    # Plot element numbers
    if elementnumbers:
        for elem in range(mesh.nelem):
            plt.text(np.mean(x[:, elem]), np.mean(y[:, elem]), str(elem + 1),
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
    if 'curves' in kwargs:
        curves = kwargs['curves']
    if 'nodemarks' in kwargs:
        nodemarks = kwargs['nodemarks']
    if 'nodenumbers' in kwargs:
        nodenumbers = kwargs['nodenumbers']
    if 'elementnumbers' in kwargs:
        elementnumbers = kwargs['elementnumbers']
    if 'curvenumbers' in kwargs:
        curvenumbers = kwargs['curvenumbers']
    if 'pointnumbers' in kwargs:
        pointnumbers = kwargs['pointnumbers']

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
                     'C' + str(curve + 1), fontsize=8, color='r',
                     ha='center', va='center')

    # Plot point numbers
    if pointnumbers:
        for point in range(mesh.npoints):
            node = mesh.points[point]
            plt.text(mesh.coor[node, 0], mesh.coor[node, 1],
                     'P' + str(point + 1), fontsize=8, color='b')

    plt.show()


def plot_over_line(mesh_pv, problem, u, points, physq=0, degfd=0, npoints=200,
                   plot_mesh=False):
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

    Keyword arguments
    -----------------
    kwargs : dict, optional
        Additional keyword arguments to pass to the plotter.add_mesh function.

    """

    # Optional arguments
    physq = kwargs.get('physq', 0)
    window_size = kwargs.get('window_size', (800, 400))

    kwargs.pop('physq', None)
    kwargs.pop('degfd', None)
    kwargs.pop('window_size', None)

    ndf = problem.elementdof[0, physq]

    assert (problem.elementdof[0, physq] == 2)
    assert (np.all(problem.elementdof[0, physq] == ndf))

    mesh_pv_plot = fill_mesh_pv(mesh_pv, problem, u, physq, degfd=[0, 1])

    glyphs = mesh_pv_plot.glyph(orient="u", scale=True, factor=0.1)

    plotter = pv.Plotter(window_size=window_size)
    # plotter.add_mesh(glyphs, show_scalar_bar=False, lighting=False,
    #                  cmap='coolwarm')
    plotter.add_mesh(mesh_pv_plot, color="lightgrey")
    plotter.add_mesh(glyphs, color="black")
    plotter.camera_position = 'xy'
    plotter.add_text((f'sol physq = {physq:d}  degfd = {physq:d}'),
                     font_size=12)
    plotter.show()
