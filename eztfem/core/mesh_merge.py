import numpy as np
from .class_mesh import Mesh


def mesh_merge(mesh1, mesh2, **kwargs):
    """
    Create a new mesh by merging two other meshes.

    Parameters
    ----------
    mesh1 : Mesh
        The first mesh to be merged.
    mesh2 : Mesh
        The second mesh to be merged.

    Keyword arguments
    -----------------
    point1, point2 : array_like
        Vectors of equal size giving the points in both meshes that need
        to be merged into one.
        NOTE: points2[i] corresponds to points1[i].

    curves1, curves2 : array_like
        Vectors of equal size giving the curves in both meshes that need
        to be merged into one.
        NOTE: curves2[i] corresponds to curves1[i].

    dir_curves2 : array_like
        Array with value of 1 or -1, indicating the direction in which
        curves2 is traversed. Must be present if curves2 is present, and of
        same length as curves2
        The curve number is curves2[i] and the curve is traversed in the
        opposite direction if dir_curves2[i] < 0.

    deletepoints1, deletepoints2 : array_like
        Vectors to indicate which points of mesh1 and mesh2 must not be
        included in the new merged mesh.
        NOTE: all double points are removed automatically already.

    deletecurves1, deletecurves2 : array_like
        Vectors to indicate which curves of mesh1 and mesh2 must not be
        included in the new merged mesh.
        NOTE: merged curves are removed automatically already.

    Returns
    -------
    mesh : Mesh
        The new merged mesh structure.

    Examples
    --------
    >>> mesh = mesh_merge(mesh1, mesh2, curves1=[1, 2], curves2=[2, -3])
    Creates a mesh from two meshes where the curves 1 and 2 of mesh1 are
    merged with curves 2 and 3 of mesh2 into two curves of the new mesh.
    Curve 3 of mesh2 is traversed in the opposite direction.

    """

    if mesh1.ndim != mesh2.ndim:
        raise ValueError('ndim different in mesh1 and mesh2')

    if mesh1.elshape != mesh2.elshape:
        raise ValueError('elshape different in mesh1 and mesh2')

    points1 = kwargs.get('points1', None)
    points2 = kwargs.get('points2', None)
    curves1 = kwargs.get('curves1', None)
    curves2 = kwargs.get('curves2', None)
    dir_curves2 = kwargs.get('dir_curves2', None)
    deletepoints1 = kwargs.get('deletepoints1', None)
    deletepoints2 = kwargs.get('deletepoints2', None)
    deletecurves1 = kwargs.get('deletecurves1', None)
    deletecurves2 = kwargs.get('deletecurves2', None)

    # booleans to indicate if the optional argument was present
    pnts1_present = points1 is not None
    pnts2_present = points2 is not None
    crvs1_present = curves1 is not None
    crvs2_present = curves2 is not None
    dir_crvs2_present = dir_curves2 is not None
    delpnts1_present = deletepoints1 is not None
    delpnts2_present = deletepoints2 is not None
    delcrvs1_present = deletecurves1 is not None
    delcrvs2_present = deletecurves2 is not None

    for option in kwargs:
        if option not in ['points1', 'points2', 'curves1', 'curves2',
                          'dir_curves2', 'deletepoints1', 'deletepoints2',
                          'deletecurves1', 'deletecurves2']:
            raise ValueError(f'Invalid option: {option}')

    if pnts1_present and pnts2_present:
        if len(points1) != len(points2):
            raise ValueError('points1 and points2 need to be of the same \
                             length')
        pnts = 1
    elif pnts1_present or pnts2_present:
        raise ValueError('both points1 and points2 need to be present')
    else:
        pnts = 0

    if crvs1_present and crvs2_present:
        if len(curves1) != len(curves2):
            raise ValueError('curves1 and curves2 need to be of the same \
                             length')
        for i in range(len(curves1)):
            if (mesh1.curves[curves1[i]].nnodes
               != mesh2.curves[abs(curves2[i])].nnodes):
                raise ValueError('number of nodes different for curves1 and \
                                 curves2')
        crvs = 1
    elif crvs1_present or crvs2_present:
        raise ValueError('both curves1 and curves2 need to be present')
    else:
        crvs = 0

    if crvs2_present and not dir_crvs2_present:
        raise ValueError('curves2 and dir_curves2 must be both present')
    elif crvs2_present and dir_crvs2_present:
        if len(curves2) != len(dir_curves2):
            raise ValueError('curves2 and dir_curves2 need to be of the same \
                             length')
        for ii, val in enumerate(curves2):
            if val < 0:
                raise ValueError('curves2 must be a non-negative curve number')
            if dir_curves2[ii] not in [1, -1]:
                raise ValueError('curves2_dir must be 1 or -1')

    # work storage for new numbering of new nodes for mesh2
    # NOTE: -1 indicates that a node has not been assigned yet
    work = np.zeros(mesh2.nnodes, dtype=int) - 1

    # arrays for deleting points and curves when merging
    delpoints1 = np.zeros(mesh1.npoints, dtype=int)
    delpoints2 = np.zeros(mesh2.npoints, dtype=int)
    delcurves1 = np.zeros(mesh1.ncurves, dtype=int)
    delcurves2 = np.zeros(mesh2.ncurves, dtype=int)

    # points
    if pnts:
        delpoints2[points2] = 1  # always delete points in points2

        # node in points2 can be removed == node in points1
        work[mesh2.points[points2]] = mesh1.points[points1]

    # delete points
    if delpnts1_present:
        delpoints1[deletepoints1] = 1

    if delpnts2_present:
        delpoints2[deletepoints2] = 1

    numdelp1 = np.sum(delpoints1)
    numdelp2 = np.sum(delpoints2)

    # curves
    if crvs:
        delcurves2[np.abs(curves2)] = 1  # always delete curves in curves2

        # nodes on curve2 can be removed == nodes on curve1
        for crv in range(len(curves2)):
            crv1 = curves1[crv]
            crv2 = curves2[crv]
            dir2 = dir_curves2[crv]
            if dir2 > 0:
                work[mesh2.curves[crv2].nodes] = mesh1.curves[crv1].nodes
            else:
                nnodes_crv2 = mesh2.curves[crv2].nnodes
                work[mesh2.curves[crv2].nodes[nnodes_crv2 - 1::-1]] \
                    = mesh1.curves[crv1].nodes

    # delete curves
    if delcrvs1_present:
        delcurves1[deletecurves1] = 1

    if delcrvs2_present:
        delcurves2[deletecurves2] = 1

    numdelc1 = np.sum(delcurves1)
    numdelc2 = np.sum(delcurves2)

    nn = mesh1.nnodes

    # new node numbers for mesh2
    for nodenr in range(mesh2.nnodes):
        if work[nodenr] > -1:
            continue  # has new node number from mesh1
        work[nodenr] = nn
        nn += 1

    # fill new mesh
    mesh = Mesh(ndim=mesh1.ndim, nnodes=nn, elshape=mesh1.elshape,
                nelem=mesh1.nelem+mesh2.nelem, elnumnod=mesh1.elnumnod,
                topology=np.zeros((mesh1.elnumnod, mesh1.nelem+mesh2.nelem),
                                  dtype=int),
                coor=np.zeros((nn, mesh1.ndim)))

    # copy elements from mesh1
    mesh.topology[:, :mesh1.nelem] = mesh1.topology

    # copy elements from mesh2
    mesh.topology[:, mesh1.nelem:] = work[mesh2.topology]

    # fill points

    # remove all double points
    numdelp2_added = 0
    for pnt2 in range(mesh2.npoints):
        if delpoints2[pnt2] == 1:
            continue  # already deleted
        if mesh2.points[pnt2] == -1:
            continue  # point not connected

        # remove double points
        if work[mesh2.points[pnt2]] in mesh1.points:
            delpoints2[pnt2] = 1
            numdelp2_added += 1

    numdelp2 += numdelp2_added
    mesh.npoints = mesh1.npoints + mesh2.npoints - numdelp1 - numdelp2
    sp = mesh1.npoints - numdelp1

    mesh.points = np.zeros(max(len(mesh1.points), mesh.npoints), dtype=int)

    # copy points mesh1
    pnt = 0
    for point in range(mesh1.npoints):
        if delpoints1[point] == 1:
            continue  # no copy of point

        mesh.points[pnt] = mesh1.points[point]
        pnt += 1

    # copy points mesh2
    pnt = 0
    for point in range(mesh2.npoints):
        if delpoints2[point] == 1:
            continue  # no copy of point
        # use new node numbers
        mesh.points[sp + pnt] = work[mesh2.points[point]]
        pnt += 1

    # fill curves
    mesh.ncurves = mesh1.ncurves + mesh2.ncurves - numdelc1 - numdelc2
    sp = mesh1.ncurves - numdelc1

    # copy curves mesh1
    # crv = 0

    for curve in range(mesh1.ncurves):
        if delcurves1[curve] == 1:
            continue  # no copy of curve

        mesh.curves.append(mesh1.curves[curve])
        # crv += 1

    # copy curves mesh2

    for curve in range(mesh2.ncurves):
        if delcurves2[curve] == 1:
            continue  # no copy of curve

        mesh.curves.append(mesh2.curves[curve])
        mesh.curves[-1].nodes = work[mesh2.curves[curve].nodes]
        # Assuming mesh.curves[sp+crv].topology is a 2D numpy array
        mesh.curves[-1].topology[:, :, 1] \
            = work[mesh2.curves[curve].topology[:, :, 1]]
        # crv += 1

    # coordinates

    # coordinates mesh1
    mesh.coor[:mesh1.nnodes] = mesh1.coor

    # coordinates mesh2
    for i in range(mesh2.nnodes):
        if work[i] < mesh1.nnodes:
            continue
        mesh.coor[work[i]] = mesh2.coor[i]

    return mesh
