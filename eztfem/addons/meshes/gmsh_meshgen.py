"""
Module for Gmsh mesh generation.
"""
from __future__ import annotations

import typing

import numpy as np
import gmsh

from ...core.meshgen import Geometry, Mesh

_GMSH_ELEMENT_MAP = {
    "tri3": {"gmsh_type": 2, "elshape": 3, "elnumnod": 3},
    "tri6": {"gmsh_type": 9, "elshape": 4, "elnumnod": 6},
    "quad4": {"gmsh_type": 3, "elshape": 5, "elnumnod": 4},
    "quad9": {"gmsh_type": 10, "elshape": 6, "elnumnod": 9},
}
_GMSH_TYPE_TO_ELEMENT = {
    value["gmsh_type"]: key for key, value in _GMSH_ELEMENT_MAP.items()
}


def gmsh_inv_index(gmsh_type: int) -> np.ndarray | None:
    """
    Return 0-based permutation mapping Gmsh local node ordering
    to internal FEM ordering.
    """
    return {
        2:  np.array([0, 2, 1], dtype=int),          # 3-node line
        9:  np.array([0, 3, 1, 4, 2, 5], dtype=int),  # 6-node triangle
        10: np.array([0, 4, 1, 5, 2, 6, 3, 7, 8], dtype=int),  # 9-node quad
    }.get(gmsh_type)


def _gmsh_element_info(element_type: str) -> dict[str, int]:
    if element_type not in _GMSH_ELEMENT_MAP:
        valid = ", ".join(sorted(_GMSH_ELEMENT_MAP))
        raise ValueError(f"Unsupported element_type '{element_type}'."
                         f"Options: {valid}.")
    return _GMSH_ELEMENT_MAP[element_type]


def _gmsh_line_element_info(element_type: int) -> tuple[int, int]:
    name, _, _, num_nodes, _, _ = (
        gmsh.model.mesh.getElementProperties(element_type)
    )
    if num_nodes == 2:
        return 1, 2
    if num_nodes == 3:
        return 2, 3
    raise ValueError(
        "Only 2-node and 3-node line elements supported for boundary curves."
        f" Got '{name}' with {num_nodes} nodes."
    )


def _gmsh_nodes_to_mesh() -> tuple[np.ndarray, dict[int, int]]:
    node_tags, coords, _ = gmsh.model.mesh.getNodes()
    coords = np.array(coords, dtype=float).reshape(-1, 3)
    coor = coords[:, :2]
    node_tag_to_index = {int(tag): idx for idx, tag in enumerate(node_tags)}
    return coor, node_tag_to_index


def _gmsh_elements_to_topology(
    node_tag_to_index: dict[int, int],
    element_type: str,
) -> tuple[np.ndarray, int, int]:
    info = _gmsh_element_info(element_type)
    gmsh_type = info["gmsh_type"]
    elshape = info["elshape"]
    elnumnod = info["elnumnod"]

    element_types, _, element_node_tags = gmsh.model.mesh.getElements(dim=2)
    element_types = np.array(element_types, dtype=int)
    matches = np.where(element_types == gmsh_type)[0]
    if matches.size == 0:
        raise ValueError(f"Gmsh mesh does not contain element type "
                         f"'{element_type}'.")

    idx = int(matches[0])

    nodes = np.array(element_node_tags[idx], dtype=int).reshape(-1, elnumnod)

    perm = gmsh_inv_index(gmsh_type)
    if perm is not None:
        nodes = nodes[:, perm]

    topology = np.vectorize(node_tag_to_index.get)(nodes).T
    return topology, elshape, elnumnod


def _gmsh_infer_element_type_or_raise() -> str:
    """
    Infer the unique supported 2D element type present in the mesh.

    Raises if:
    - no supported 2D element types are present, or
    - more than one supported 2D element type is present.
    """
    element_types, _, _ = gmsh.model.mesh.getElements(dim=2)
    present = sorted(set(int(t) for t in element_types))
    supported_present = [t for t in present if t in _GMSH_TYPE_TO_ELEMENT]

    if not supported_present:
        raise ValueError(
            "No supported 2D elements found in the Gmsh mesh. "
            f"Present gmsh types: {present}"
        )

    if len(supported_present) == 1:
        return _GMSH_TYPE_TO_ELEMENT[supported_present[0]]

    names = ", ".join(_GMSH_TYPE_TO_ELEMENT[t] for t in supported_present)
    raise ValueError(
        "Multiple supported 2D element types found in the Gmsh mesh; "
        "this converter expects a single type. "
        f"Found: {names} (gmsh types: {supported_present})."
    )


def _gmsh_physical_point_tags() -> list[int]:
    return [int(tag) for _, tag in gmsh.model.getPhysicalGroups(dim=0)]


def _gmsh_points_from_physical_groups(
    node_tag_to_index: dict[int, int],
) -> np.ndarray:
    """
    Collect unique mesh node indices belonging to physical point group (dim=0).
    Returns an array of global node indices (into mesh.coor).
    """
    point_node_indices: set[int] = set()

    for physical_tag in _gmsh_physical_point_tags():
        entities = gmsh.model.getEntitiesForPhysicalGroup(0, physical_tag)
        for entity_tag in entities:
            node_tags, _, _ = gmsh.model.mesh.getNodes(dim=0, tag=entity_tag)
            for nt in node_tags:
                point_node_indices.add(node_tag_to_index[int(nt)])

    return np.array(sorted(point_node_indices), dtype=int)


def _gmsh_physical_curve_tags() -> list[int]:
    return [int(tag) for _, tag in gmsh.model.getPhysicalGroups(dim=1)]


def _gmsh_curve_geometry(
    node_tag_to_index: dict[int, int],
    physical_tag: int,
) -> Geometry | None:
    entities = gmsh.model.getEntitiesForPhysicalGroup(1, physical_tag)
    if not entities:
        return None

    curve_nodes: list[int] = []
    curve_node_lookup: dict[int, int] = {}
    local_topology: list[list[int]] = []
    line_elshape = None
    line_elnumnod = None

    for entity_tag in entities:
        element_types, _, element_node_tags = gmsh.model.mesh.getElements(
            dim=1, tag=entity_tag
        )
        for elem_type, node_tags in zip(element_types, element_node_tags):
            elshape, elnumnod = _gmsh_line_element_info(int(elem_type))
            if line_elnumnod is None:
                line_elshape = elshape
                line_elnumnod = elnumnod
            elif line_elnumnod != elnumnod:
                raise ValueError("Mixed line element orders not supported.")

            nodes = np.array(node_tags, dtype=int).reshape(-1, elnumnod)

            perm = gmsh_inv_index(int(elem_type))
            if perm is not None:
                nodes = nodes[:, perm]

            for element_nodes in nodes:
                local_nodes: list[int] = []
                for node_tag in element_nodes:
                    global_index = node_tag_to_index[int(node_tag)]
                    local_index = curve_node_lookup.get(global_index)
                    if local_index is None:
                        local_index = len(curve_nodes)
                        curve_node_lookup[global_index] = local_index
                        curve_nodes.append(global_index)
                    local_nodes.append(local_index)
                local_topology.append(local_nodes)

    if line_elnumnod is None or line_elshape is None:
        return None

    curve_nodes_arr = np.array(curve_nodes, dtype=int)
    local_topology_arr = np.array(local_topology, dtype=int).T

    topology = np.zeros((line_elnumnod, local_topology_arr.shape[1], 2),
                        dtype=int)
    topology[:, :, 0] = local_topology_arr
    topology[:, :, 1] = curve_nodes_arr[local_topology_arr]

    curve = Geometry(
        elshape=line_elshape,
        ndim=2,
        elnumnod=line_elnumnod,
        nnodes=curve_nodes_arr.size,
        nelem=local_topology_arr.shape[1],
        topology=topology,
        nodes=curve_nodes_arr,
    )
    return curve


def gmsh_mesh2d(
    model_builder: typing.Callable[[typing.Any], None],
    gmsh_options: dict[str, float] | None = None
) -> Mesh:
    """
    Build a 2D mesh in Gmsh and convert it to an eztfem Mesh.

    - Calls `model_builder(open_gui=False)` to define the geometry/model.
    - Generates a 2D mesh.
    - Infers the unique supported 2D element type present in the mesh.
    - Raises if multiple supported 2D element types are present.

    Parameters
    ----------
    model_builder : callable
        Function that receives the gmsh module (or compatible signature) and
        defines the model. Must accept open_gui keyword.

    Returns
    -------
    mesh : Mesh
        The generated mesh in eztfem format.
    """
    should_finalize = False
    if not gmsh.isInitialized():
        gmsh.initialize()
        should_finalize = True

    if gmsh_options:
        for k, v in gmsh_options.items():
            gmsh.option.setNumber(k, v)

    try:
        model_builder(open_gui=False)
        gmsh.model.mesh.generate(dim=2)

        resolved_type = _gmsh_infer_element_type_or_raise()

        coor, node_tag_to_index = _gmsh_nodes_to_mesh()
        topology, elshape, elnumnod = _gmsh_elements_to_topology(
            node_tag_to_index, resolved_type
        )

        points = _gmsh_points_from_physical_groups(node_tag_to_index)

        mesh = Mesh(
            ndim=2,
            nnodes=coor.shape[0],
            nelem=topology.shape[1],
            elshape=elshape,
            elnumnod=elnumnod,
            topology=topology,
            coor=coor,
            npoints=int(points.size),
            points=points,
            curves=[],
        )

        curves: list[Geometry] = []
        for physical_tag in _gmsh_physical_curve_tags():
            curve = _gmsh_curve_geometry(node_tag_to_index, physical_tag)
            if curve is not None:
                curves.append(curve)

        mesh.curves = curves
        mesh.ncurves = len(curves)

        return mesh
    finally:
        if should_finalize:
            gmsh.finalize()
