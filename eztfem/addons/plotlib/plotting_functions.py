"""Module for plotting functions using PyVista."""
import copy
import typing

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import pyvista as pv
from pyvista import CellType

from ...core.gauss import gauss_legendre
from ...core.meshgen import Mesh
from ...core.pos_array import pos_array, pos_array_vec
from ...core.problem import Problem
from ...core.shapefunc import basis_function
from ...core.vector import Vector

FloatArray: typing.TypeAlias = npt.NDArray[np.floating]
IntArray: typing.TypeAlias = npt.NDArray[np.integer]


def fill_mesh_pv(
    mesh_pv: pv.PolyData,
    problem: Problem,
    u: FloatArray,
    physq: int,
    degfd: int | typing.Sequence[int],
) -> pv.PolyData:
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


def fill_mesh_pv_vector(
    mesh_pv: pv.PolyData,
    problem: Problem,
    vector: Vector,
    degfd: int,
) -> pv.PolyData:
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


def plot_mesh_pv(
    mesh_pv: pv.PolyData,
    *,
    style: typing.Literal[
        'surface', 'wireframe', 'points', 'points_gaussian'
    ] = "wireframe",
    color: str = "black",
    window_size: tuple[int, int] = (800, 400),
    **kwargs: typing.Any,
) -> None:
    """
    Plot a mesh.

    Parameters
    ----------
    mesh_pv : pyvista.PolyData
        The mesh to be plotted.
    style : {'surface', 'wireframe', 'points', 'points_gaussian'}, optional
        Edge style passed to plotter.add_mesh. Default "wireframe".
    color : str, optional
        Edge color passed to plotter.add_mesh. Default "black".
    window_size : tuple of int, optional
        Size of the PyVista rendering window. Default (800, 400).
    kwargs : dict, optional
        Additional keyword arguments to pass to the plotter.add_mesh function.

    Notes
    -----
    This tutorial was followed to hide the internal edges:
    https://github.com/pyvista/pyvista/discussions/5777

    """
    surface = mesh_pv.separate_cells().extract_surface(
        nonlinear_subdivision=4, algorithm='dataset_surface')
    edges = surface.extract_feature_edges()

    plotter = pv.Plotter(window_size=list(window_size))
    plotter.add_mesh(surface)
    plotter.add_mesh(edges, style=style, color=color, **kwargs)
    plotter.camera_position = 'xy'
    plotter.show()


def plot_sol(
    mesh_pv: pv.PolyData,
    problem: Problem,
    u: FloatArray,
    *,
    physq: int = 0,
    degfd: int = 0,
    window_size: tuple[int, int] = (800, 400),
    **kwargs: typing.Any,
) -> None:
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
    physq : int, optional
        Physical quantity number. Default 0.
    degfd : int, optional
        Degree of freedom to plot. Default 0.
    window_size : tuple of int, optional
        Size of the PyVista rendering window. Default (800, 400).
    kwargs : dict, optional
        Additional keyword arguments to pass to the plotter.add_mesh function.

    """
    mesh_pv_plot = fill_mesh_pv(mesh_pv, problem, u, physq, degfd)

    plotter = pv.Plotter(window_size=list(window_size))
    plotter.add_mesh(mesh_pv_plot, scalars="u", **kwargs)
    plotter.camera_position = 'xy'
    plotter.add_text((f'sol physq = {physq:d}  degfd = {physq:d}'),
                     font_size=12)
    plotter.show()


def plot_sol_contour(
    mesh_pv: pv.PolyData,
    problem: Problem,
    u: FloatArray,
    nlevels: int = 10,
    *,
    physq: int = 0,
    degfd: int = 0,
    window_size: tuple[int, int] = (800, 400),
    **kwargs: typing.Any,
) -> None:
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
    physq : int, optional
        Physical quantity number. Default 0.
    degfd : int, optional
        Degree of freedom to plot. Default 0.
    window_size : tuple of int, optional
        Size of the PyVista rendering window. Default (800, 400).
    kwargs : dict, optional
        Additional keyword arguments to pass to the plotter.add_mesh function.

    """
    mesh_pv_plot = fill_mesh_pv(mesh_pv, problem, u, physq, degfd)

    contours = mesh_pv_plot.contour(nlevels, scalars='u')

    plotter = pv.Plotter(window_size=list(window_size))
    plotter.add_mesh(mesh_pv_plot, color="lightgrey", **kwargs)
    plotter.add_mesh(contours)  # color="black", line_width=1)
    plotter.camera_position = 'xy'
    plotter.add_text((f'sol physq = {physq:d}  degfd = {physq:d}'),
                     font_size=12)
    plotter.show()


def plot_vector(
    mesh_pv: pv.PolyData,
    problem: Problem,
    vector: Vector,
    degfd: int = 0,
    *,
    window_size: tuple[int, int] = (800, 400),
    **kwargs: typing.Any,
) -> None:
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
    window_size : tuple of int, optional
        Size of the PyVista rendering window. Default (800, 400).
    kwargs : dict, optional
        Additional keyword arguments to pass to the plotter.add_mesh function.

    """
    mesh_pv_plot = fill_mesh_pv_vector(mesh_pv, problem, vector, degfd)

    plotter = pv.Plotter(window_size=list(window_size))
    plotter.add_mesh(mesh_pv_plot, scalars="u", **kwargs)
    plotter.camera_position = 'xy'
    plotter.add_text((f'sol vec = {vector.vec:d}  degfd = {degfd:d}'),
                     font_size=12)
    plotter.show()


def plot_vector_contours(
    mesh_pv: pv.PolyData,
    problem: Problem,
    vector: Vector,
    degfd: int = 0,
    nlevels: int = 10,
    *,
    window_size: tuple[int, int] = (800, 400),
    **kwargs: typing.Any,
) -> None:
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
    window_size : tuple of int, optional
        Size of the PyVista rendering window. Default (800, 400).
    kwargs : dict, optional
        Additional keyword arguments to pass to the plotter.add_mesh function.

    """
    mesh_pv_plot = fill_mesh_pv_vector(mesh_pv, problem, vector, degfd)

    contours = mesh_pv_plot.contour(nlevels, scalars='u')

    plotter = pv.Plotter(window_size=list(window_size))
    plotter.add_mesh(mesh_pv_plot, color="lightgrey", **kwargs)
    plotter.add_mesh(contours)  # color="black", line_width=1)
    plotter.camera_position = 'xy'
    plotter.add_text((f'sol vec = {vector.vec:d}  degfd = {degfd:d}'),
                     font_size=12)
    plotter.show()


def plot_mesh_plt(
    mesh: Mesh,
    *,
    nodemarks: bool = False,
    nodenumbers: bool = False,
    elementnumbers: bool = False,
) -> None:
    """
    Plot mesh structure.

    Parameters
    ----------
    mesh : Mesh
        Mesh structure.
    nodemarks : bool, optional
        Plot node marks. Default False.
    nodenumbers : bool, optional
        Plot node numbers. Default False.
    elementnumbers : bool, optional
        Plot element numbers. Default False.
    """

    if mesh.ndim != 2:
        raise ValueError('Only 2D meshes can be plotted')

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


def plot_curves(
    mesh: Mesh,
    *,
    curves: typing.Sequence[int] | None = None,
    nodemarks: bool = False,
    nodenumbers: int = 0,
    elementnumbers: bool = False,
    curvenumbers: bool = False,
    pointnumbers: bool = False,
) -> None:
    """
    Plot curves (and optionally points) of the mesh.

    Parameters
    ----------
    mesh : Mesh
        Mesh structure.
    curves : sequence of int, optional
        The curves to be plotted. Default is all curves.
    nodemarks : bool, optional
        Plot node marks. Default False.
    nodenumbers : int, optional
        Plot node numbers (0=off, 1=local, 2=global). Default 0.
    elementnumbers : bool, optional
        Plot element numbers. Default False.
    curvenumbers : bool, optional
        Plot curve numbers. Default False.
    pointnumbers : bool, optional
        Plot point numbers. Default False.
    """

    if mesh.ndim != 2:
        raise ValueError('Only 2D curves can be plotted')

    if curves is None:
        curves = list(range(mesh.ncurves))

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


def plot_points_curves(mesh: Mesh) -> None:
    """Plot all points and curves .

    Parameters
    ----------
    mesh : Mesh
        Mesh object containing the curves and points to be plotted.
    """

    plot_curves(mesh, pointnumbers=True, curvenumbers=True)


def plot_sol_over_line(
    mesh_pv: pv.PolyData,
    problem: Problem,
    u: FloatArray,
    points: typing.Sequence[FloatArray],
    physq: int = 0,
    degfd: int = 0,
    npoints: int = 200,
    plot_mesh: bool = False,
) -> pv.PolyData:
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

    # PyVista's line/sampling helpers expect float64 arrays specifically.
    point0 = np.asarray(points[0], dtype=np.float64)
    point1 = np.asarray(points[1], dtype=np.float64)

    if plot_mesh:
        line = pv.Line(point0, point1)
        p = pv.Plotter(window_size=[800, 400])
        p.add_mesh(mesh_pv, color="w")
        p.add_mesh(line, color="b")
        p.camera_position = 'xy'
        p.show()

    mesh_pv_plot = fill_mesh_pv(mesh_pv, problem, u, physq, degfd)

    mesh_pv_plot.plot_over_line(point0, point1, resolution=npoints)

    return mesh_pv_plot.sample_over_line(point0, point1,
                                         resolution=npoints)


def plot_vector_over_line(
    mesh_pv: pv.PolyData,
    problem: Problem,
    vector: Vector,
    points: typing.Sequence[FloatArray],
    degfd: int = 0,
    npoints: int = 200,
    plot_mesh: bool = False,
) -> pv.PolyData:
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

    # PyVista's line/sampling helpers expect float64 arrays specifically.
    point0 = np.asarray(points[0], dtype=np.float64)
    point1 = np.asarray(points[1], dtype=np.float64)

    if plot_mesh:
        line = pv.Line(point0, point1)
        p = pv.Plotter(window_size=[800, 400])
        p.add_mesh(mesh_pv, color="w")
        p.add_mesh(line, color="b")
        p.camera_position = 'xy'
        p.show()

    mesh_pv_plot = fill_mesh_pv_vector(mesh_pv, problem, vector, degfd)

    mesh_pv_plot.plot_over_line(point0, point1, resolution=npoints)

    return mesh_pv_plot.sample_over_line(point0, point1,
                                         resolution=npoints)


def plot_quiver(
    mesh_pv: pv.PolyData,
    problem: Problem,
    u: FloatArray,
    *,
    physq: int = 0,
    window_size: tuple[int, int] = (800, 400),
    scale: float = 0.1,
    **kwargs: typing.Any,
) -> None:
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
    physq : int, optional
        Physical quantity number. Default 0.
    window_size : tuple of int, optional
        Size of the PyVista rendering window. Default (800, 400).
    scale : float, optional
        Scaling factor for the arrows (default = 0.1)
    kwargs : dict, optional
        Additional keyword arguments to pass to the plotter.add_mesh function.

    """
    ndf = problem.elementdof[0, physq]

    assert (problem.elementdof[0, physq] == 2)
    assert np.all(problem.elementdof[0, physq] == ndf)

    mesh_pv_plot = fill_mesh_pv(mesh_pv, problem, u, physq, degfd=[0, 1])

    glyphs = mesh_pv_plot.glyph(orient="u", scale=True, factor=scale)

    plotter = pv.Plotter(window_size=list(window_size))
    # plotter.add_mesh(glyphs, show_scalar_bar=False, lighting=False,
    #                  cmap='coolwarm')
    plotter.add_mesh(mesh_pv_plot, color="lightgrey")
    plotter.add_mesh(glyphs, color="black")
    plotter.camera_position = 'xy'
    plotter.add_text((f'sol physq = {physq:d}  degfd = {physq:d}'),
                     font_size=12)
    plotter.show()


def plot_basis_function(
    shape: typing.Literal['quad', 'triangle'],
    intpol: typing.Literal[
        'P0', 'P1', 'P1+', 'P2', 'P2+', 'Q1', 'Q1+', 'Q2'
    ],
    degfd: int,
    *,
    n: int = 10,
    plot3d: bool = True,
    edges: bool = True,
    axisoff: bool = False,
    show: bool = True,
    window_size: tuple[int, int] = (800, 400),
    **kwargs: typing.Any,
) -> tuple[pv.Plotter, pv.UnstructuredGrid]:
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

    plotter = pv.Plotter(window_size=list(window_size))

    mesh_kwargs = {
        'scalars': 'phi',
        'show_edges': edges,
    }
    mesh_kwargs.update(kwargs)

    # add_mesh has many precisely-typed keyword arguments that a generic
    # forwarded-kwargs dict can't satisfy statically.
    plotter.add_mesh(mesh, **mesh_kwargs)  # type: ignore[arg-type]

    if not plot3d:
        plotter.view_xy()  # type: ignore[call-arg]

    title = f"phi      node = {degfd_idx + 1}   intpol = {intpol}"
    plotter.add_title(title)

    if axisoff:
        plotter.remove_scalar_bar()  # type: ignore[call-arg]
    else:
        plotter.show_bounds()  # type: ignore[call-arg]

    if show:
        plotter.show()

    return plotter, mesh


def _reference_quad_mesh(
    n: int,
) -> tuple[FloatArray, IntArray, npt.NDArray[np.uint8]]:
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


def _reference_triangle_mesh(
    n: int,
) -> tuple[FloatArray, IntArray, npt.NDArray[np.uint8]]:
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
    celltype_list = []

    np_idx = 0
    for j in range(n):
        first_count = n - j
        for i in range(first_count):
            n0 = np_idx + i
            n1 = np_idx + i + 1
            n2 = np_idx + i + n - j + 1
            cells.extend([3, n0, n1, n2])
            celltype_list.append(CellType.TRIANGLE)

        second_count = n - j - 1
        for i in range(second_count):
            n0 = np_idx + i + 1
            n1 = np_idx + i + n - j + 2
            n2 = np_idx + i + n - j + 1
            cells.extend([3, n0, n1, n2])
            celltype_list.append(CellType.TRIANGLE)

        np_idx += n - j + 1

    cell_array = np.array(cells, dtype=np.int64)
    celltypes = np.array(celltype_list, dtype=np.uint8)

    return points, cell_array, celltypes


def plot_gauss_legendre(
    shape: str,
    *,
    n: int | None = None,
    p: int | None = None,
    ax: plt.Axes | None = None,
    marker: str = '+',
    color: str = 'k',
    markersize: float = 8,
    show: bool = True,
    **kwargs: typing.Any,
) -> plt.Axes:
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
