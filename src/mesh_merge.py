from src.mesh_class import Mesh
import numpy as np

def mesh_merge(mesh1, mesh2, **kwargs):

    # Test meshes
    if mesh1.ndim != mesh2.ndim:
        raise ValueError('ndim different in mesh1 and mesh2')
    if mesh1.elshape != mesh2.elshape:
        raise ValueError('elshape different in mesh1 and mesh2')

    def check_type(var):
        var_out = [var] if isinstance(var, (int, np.integer)) else var
        return var_out

    # Default values for optional arguments
    points1 = check_type(kwargs.get('points1', []))
    points2 = check_type(kwargs.get('points2', []))
    curves1 = check_type(kwargs.get('curves1', []))
    curves2 = check_type(kwargs.get('curves2', []))
    deletepoints1 = check_type(kwargs.get('deletepoints1', []))
    deletepoints2 = check_type(kwargs.get('deletepoints2', []))
    deletecurves1 = check_type(kwargs.get('deletecurves1', []))
    deletecurves2 = check_type(kwargs.get('deletecurves2', []))

    # if len(points1) != len(points2) and len(points1) != 0:
    #     raise ValueError('points1 and points2 need to be of the same length or one should be empty')
    
    # # Check curves
    # if len(curves1) != len(curves2) and len(curves1) != 0:
    #     raise ValueError('curves1 and curves2 need to be of the same length or one should be empty')
    
    # Handle points
    if points1 is None:  # If no points provided, we assume all points are to be merged
        points1 = list(range(len(mesh1.points)))
        points2 = list(range(len(mesh2.points)))

    # Handle deletepoints
    if deletepoints1 is None:
      points1 = list(set(points1))
    else:
      points1 = list(set(points1) - set(deletepoints1))

    if deletepoints2 is None:
      points2 = list(set(points2))
    else:
      points2 = list(set(points2) - set(deletepoints2))

    merged_mesh = Mesh()
    
    merged_mesh.points = [mesh1.points[i] for i in points1] + [mesh2.points[i] for i in points2]
    
    # Handle curves
    if curves1 is None:  # If no curves provided, we assume all curves are to be merged
        curves1 = list(range(len(mesh1.curves)))
        curves2 = list(range(len(mesh2.curves)))

    # Handle deletecurves
    if deletecurves1 is None:
      curves1 = list(set(curves1) )
    else:
      curves1 = list(set(curves1) - set(deletecurves1))
       
    if deletecurves2 is None:
      curves2 = list(set(curves2) )
    else:
      curves2 = list(set(curves2) - set(deletecurves2))    

    merged_mesh.curves = [mesh1.curves[i] for i in curves1] + [mesh2.curves[i] for i in curves2]

    return merged_mesh
