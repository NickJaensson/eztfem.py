import typing
import pyvista as pv
import numpy as np

if typing.TYPE_CHECKING:
    from ...core.meshgen import Mesh


# TODO: Make a self-documenting int-enum for elshape
def find_elshape_vtk(elshape: int):
    """
    Find the VTK element shape and index mapping for a given element shape.

    Parameters
    ----------
    elshape : int
        The element shape identifier (see Mesh class definition).

    Returns
    -------
    elshape_vtk : int
        The VTK element shape identifier.
    index_vtk : list of int
        The index mapping for the VTK element ordering.

    Raises
    ------
    ValueError
        If the element shape identifier is invalid.
    """

    match elshape:
        case 1:  # 2-node line elements.
            elshape_vtk = 3
            index_vtk = [0, 1]
        case 2:  # 3-node line elements.
            elshape_vtk = 21
            index_vtk = [0, 2, 1]
        case 3:  # 3-node triangle.
            elshape_vtk = 5
            index_vtk = [0, 1, 2]
        case 4:  # 6-node triangle.
            elshape_vtk = 22
            index_vtk = [0, 2, 4, 1, 3, 5]
        case 5:  # 4-node quadrilateral.
            elshape_vtk = 9
            index_vtk = [0, 1, 2, 3]
        case 6:  # 9-node quadrilateral.
            elshape_vtk = 28
            index_vtk = [0, 2, 4, 6, 1, 3, 5, 7, 8]
        case 7:  # 7-node triangle.
            elshape_vtk = 34
            index_vtk = [0, 2, 4, 1, 3, 5, 6]
        case 9:  # 5-node quadrilateral.
            elshape_vtk = 9
            index_vtk = [0, 1, 2, 3, 4]
        case 10:  # 4-node triangle.
            elshape_vtk = 5
            index_vtk = [0, 1, 2, 3]
        case _:
            raise ValueError(f"invalid element shape, elshape = {elshape}")

    return elshape_vtk, index_vtk


def generate_pyvista_mesh(mesh: "Mesh"):
    """
    Generate a PyVista mesh from a given mesh object.

    Parameters
    ----------
    mesh : Mesh
        The Mesh object containing coordinates, topology, element shape, and
        other properties.

    Returns
    -------
    pv.UnstructuredGrid
        The generated PyVista unstructured grid.
    """

    # Append column of zeros to convert to 3D mesh
    coor = np.hstack((mesh.coor, np.zeros((mesh.coor.shape[0], 1))))

    # Define the cell types
    celltype, index_vtk = find_elshape_vtk(mesh.elshape)

    # Prepend the element type
    tmp1 = np.full((mesh.topology.shape[1], 1), mesh.elnumnod)
    tmp2 = mesh.topology[index_vtk, :]
    topo = np.hstack((tmp1, tmp2.T))

    celltypes = np.full(mesh.nelem, celltype)

    # Create the unstructured grid and return
    return pv.UnstructuredGrid(topo, celltypes, coor)
