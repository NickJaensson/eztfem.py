# import gmsh

def main(gmsh):

    # gmsh.initialize()
 
    gmsh.model.add("rect_with_hole")

    width  = 2.0          # rectangle width
    height = 1.2          # rectangle height

    # Create rectangle surface
    gmsh.model.occ.addRectangle(0.0, 0.0, 0.0, width, height)

    # Synchronize CAD kernel with the model
    gmsh.model.occ.synchronize()

    # Open the GUI when running this example with Python
    # gmsh.fltk.run()

    # gmsh.finalize()

if __name__ == "__main__":
    main()
