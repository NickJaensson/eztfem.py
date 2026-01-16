'''
Module for Gmsh mesh generation.
'''
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


def _gmsh_element_info(element_type: str) -> dict[str, int]:
    if element_type not in _GMSH_ELEMENT_MAP:
        valid = ", ".join(sorted(_GMSH_ELEMENT_MAP))
        raise ValueError(f"Unsupported element_type '{element_type}'. Options: {valid}.")
    return _GMSH_ELEMENT_MAP[element_type]


def _gmsh_line_element_info(element_type: int) -> tuple[int, int]:
    name, _, _, num_nodes, _, _ = gmsh.model.mesh.getElementProperties(
        element_type
    )
    if num_nodes == 2:
        return 1, 2
    if num_nodes == 3:
        return 2, 3
    raise ValueError(
        "Only 2-node and 3-node line elements are supported for boundary curves."
        f" Got '{name}' with {num_nodes} nodes."
    )


def _gmsh_nodes_to_mesh() -> tuple[np.ndarray, dict[int, int]]:
    node_tags, coords, _ = gmsh.model.mesh.getNodes()
    coords = np.array(coords, dtype=float).reshape(-1, 3)
    coor = coords[:, :2]
    node_tag_to_index = {tag: idx for idx, tag in enumerate(node_tags)}
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
        raise ValueError(
            f"Gmsh mesh does not contain element type '{element_type}'."
        )
    idx = int(matches[0])

    nodes = np.array(element_node_tags[idx], dtype=int).reshape(-1, elnumnod)
    topology = np.vectorize(node_tag_to_index.get)(nodes).T
    return topology, elshape, elnumnod


def _gmsh_detect_element_type(element_type: str | None) -> str:
    element_types, _, _ = gmsh.model.mesh.getElements(dim=2)
    element_types = np.array(element_types, dtype=int)
    supported = [gmsh_type for gmsh_type in element_types
                 if gmsh_type in _GMSH_TYPE_TO_ELEMENT]

    if element_type is not None:
        gmsh_type = _gmsh_element_info(element_type)["gmsh_type"]
        if gmsh_type not in element_types:
            raise ValueError(
                f"Gmsh mesh does not contain element type '{element_type}'."
            )
        return element_type

    unique_types = sorted(set(supported))
    if len(unique_types) == 1:
        return _GMSH_TYPE_TO_ELEMENT[unique_types[0]]
    if not unique_types:
        raise ValueError("No supported 2D elements found in the Gmsh mesh.")
    available = ", ".join(_GMSH_TYPE_TO_ELEMENT[gmsh_type]
                           for gmsh_type in unique_types)
    raise ValueError(
        "Multiple supported 2D element types found; specify element_type. "
        f"Found: {available}."
    )


def _gmsh_physical_curve_tags() -> list[int]:
    return [tag for _, tag in gmsh.model.getPhysicalGroups(dim=1)]


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
            elshape, elnumnod = _gmsh_line_element_info(elem_type)
            if line_elnumnod is None:
                line_elshape = elshape
                line_elnumnod = elnumnod
            elif line_elnumnod != elnumnod:
                raise ValueError("Mixed line element orders are not supported.")

            nodes = np.array(node_tags, dtype=int).reshape(-1, elnumnod)
            for element_nodes in nodes:
                local_nodes: list[int] = []
                for node_tag in element_nodes:
                    global_index = node_tag_to_index[node_tag]
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
    topology = np.zeros((line_elnumnod, local_topology_arr.shape[1], 2), dtype=int)
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
    *,
    element_type: str | None = None,
    msh_file: str | None = None,
    model_builder: typing.Callable[[typing.Any], None] | None = None,
    generate: bool = True,
) -> Mesh:
    """
    Load or build a 2D mesh in Gmsh and convert it to an eztfem Mesh.

    Parameters
    ----------
    element_type : {'tri3', 'tri6', 'quad4', 'quad9'}, optional
        Requested 2D element type. If None, the type is inferred when possible.
    msh_file : str, optional
        Path to a .msh file to load.
    model_builder : callable, optional
        Function that receives the gmsh module and defines the model.
    generate : bool, optional
        Whether to run gmsh.model.mesh.generate(2) after building the model.

    Returns
    -------
    mesh : Mesh
        The generated mesh in eztfem format.
    """
    if (msh_file is None) == (model_builder is None):
        raise ValueError("Specify exactly one of msh_file or model_builder.")

    should_finalize = False
    if not gmsh.isInitialized():
        gmsh.initialize()
        should_finalize = True

    try:
        if msh_file is not None:
            gmsh.open(msh_file)
        else:
            model_builder(open_gui=False)
            if generate:
                gmsh.model.mesh.generate(2)

        resolved_type = _gmsh_detect_element_type(element_type)
        coor, node_tag_to_index = _gmsh_nodes_to_mesh()
        topology, elshape, elnumnod = _gmsh_elements_to_topology(
            node_tag_to_index, resolved_type
        )

        mesh = Mesh(
            ndim=2,
            nnodes=coor.shape[0],
            nelem=topology.shape[1],
            elshape=elshape,
            elnumnod=elnumnod,
            topology=topology,
            coor=coor,
            npoints=0,
            points=np.array([], dtype=int),
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
