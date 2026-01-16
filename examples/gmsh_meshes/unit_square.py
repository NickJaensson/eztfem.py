import gmsh

def build_unit_square(open_gui: bool = True) -> None:

    if not gmsh.isInitialized():
        gmsh.initialize()

    width, height = 1.0, 1.0  # rectangle width and height

    # Create rectangle surface
    gmsh.model.occ.addRectangle(0.0, 0.0, 0.0, width, height)

    # Synchronize CAD kernel with the model
    gmsh.model.occ.synchronize()

    # Add physical points and curves
    for num in [1, 2, 3, 4]:
        gmsh.model.addPhysicalGroup(0, [num])  # add points
        gmsh.model.addPhysicalGroup(1, [num])  # add curves

    # Open the GUI when running this example with Python
    if open_gui:
        gmsh.fltk.run()

if __name__ == "__main__":
    build_unit_square()
