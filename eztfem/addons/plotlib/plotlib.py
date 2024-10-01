import numpy as np
from ...core.pos_array import pos_array
import pyvista as pv


def plot_mesh(mesh_pv, **kwargs):
    """
    Plot a mesh.

    Parameters
    ----------
    mesh_pv : pyvista.PolyData
        The mesh to be plotted.
    **kwargs : dict, optional
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
    **kwargs : dict, optional
        Additional keyword arguments to pass to the plotter.add_mesh function.
    """

    # Optional arguments
    physq = kwargs.get('physq', 0)
    degfd = kwargs.get('degfd', 0)
    window_size = kwargs.get('window_size', (800, 400))

    kwargs.pop('window_size', None)

    nnodes = mesh_pv.number_of_points

    # Set coordinates and values
    u_plot = np.zeros(nnodes)

    for node in range(nnodes):
        posn, _ = pos_array(problem, node, physq=physq, order='DN')
        u_plot[node] = u[posn[0][degfd]]

    mesh_pv.point_data['u'] = u_plot

    plotter = pv.Plotter(window_size=window_size)
    plotter.add_mesh(mesh_pv, scalars="u", **kwargs)
    plotter.camera_position = 'xy'
    plotter.add_text((f'sol physq = {physq:d}  degfd = {physq:d}'),
                     font_size=12)
    plotter.show()
