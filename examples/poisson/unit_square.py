'''
Gmsh script to build a unit square geometry.
NOTE: This script only creates the geometry and physical groups. To generate
a mesh, use the gmsh_mesh2d function in eztfem.addons.meshes.gmsh_meshgen and
pass this script as the model_builder argument.
NOTE2: This script can also be run directly with Python to open the Gmsh GUI
and visualize the geometry.
'''
import gmsh


def build_unit_square(
    open_gui: bool = True,
    element_size: float | None = None,
) -> None:
    """
    Create a unit square geometry in Gmsh and optionally open the GUI.
    """

    if not gmsh.isInitialized():
        gmsh.initialize()

    width, height = 1.0, 1.0  # rectangle width and height

    # Create rectangle surface
    gmsh.model.occ.addRectangle(0.0, 0.0, 0.0, width, height)

    # Synchronize CAD kernel with the model
    gmsh.model.occ.synchronize()

    if element_size is not None:
        for index, value in enumerate(element_size, start=1):
            gmsh.model.mesh.setSize([(0, index)], value)

    # Add physical points and curves
    for num in [1, 2, 3, 4]:
        gmsh.model.addPhysicalGroup(0, [num])  # add points
        gmsh.model.addPhysicalGroup(1, [num])  # add curves

    # Open the GUI when running this example with Python
    if open_gui:
        gmsh.fltk.run()


if __name__ == "__main__":
    build_unit_square()
