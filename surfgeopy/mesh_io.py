"""Mesh input/output helpers."""

from pathlib import Path
from typing import Dict, Iterable, Optional, Sequence, Tuple

import numpy as np

__all__ = ["read_gmsh_mesh", "write_gmsh_mesh"]

_GMSH_TRIANGLE_TYPES = {
    2: (1, 3),
    9: (2, 6),
    21: (3, 10),
    23: (4, 15),
    25: (5, 21),
}
_TRIANGLE_TYPE_BY_NODE_COUNT = {
    node_count: element_type
    for element_type, (_, node_count) in _GMSH_TRIANGLE_TYPES.items()
}


def read_gmsh_mesh(
    mesh_path: str,
    *,
    preserve_order: bool = False,
    element_types: Optional[Sequence[int]] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Read triangular surface elements from an ASCII Gmsh ``.msh`` file.

    By default high-order triangles are linearized to their three corner nodes,
    which matches Surfgeopy's reference-triangle workflow. Set
    ``preserve_order=True`` to keep all Gmsh nodes in each triangle.
    """
    path = Path(mesh_path)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {mesh_path}")

    lines = path.read_text(encoding="utf-8").splitlines()
    version = _read_mesh_version(lines)
    allowed = set(element_types) if element_types is not None else set(_GMSH_TRIANGLE_TYPES)
    unsupported = allowed.difference(_GMSH_TRIANGLE_TYPES)
    if unsupported:
        supported = ", ".join(str(element_type) for element_type in sorted(_GMSH_TRIANGLE_TYPES))
        requested = ", ".join(str(element_type) for element_type in sorted(unsupported))
        raise ValueError(
            f"Unsupported triangular Gmsh element type(s): {requested}. "
            f"Supported element types are: {supported}"
        )

    if version.startswith("2."):
        node_data, elements = _read_gmsh22(lines, allowed)
    elif version.startswith("4."):
        node_data, elements = _read_gmsh4(lines, allowed)
    else:
        raise ValueError(f"Unsupported Gmsh mesh version: {version}")

    if not elements:
        raise ValueError("Gmsh file does not contain supported triangular surface elements")

    ordered_tags = list(node_data)
    node_to_index = {tag: index for index, tag in enumerate(ordered_tags)}
    vertices = np.array([node_data[tag] for tag in ordered_tags], dtype=float)

    faces = []
    for element_type, node_tags in elements:
        _, node_count = _GMSH_TRIANGLE_TYPES[element_type]
        if len(node_tags) != node_count:
            raise ValueError(
                f"Gmsh element type {element_type} expected {node_count} nodes, "
                f"got {len(node_tags)}"
            )
        selected_tags = node_tags if preserve_order else node_tags[:3]
        try:
            faces.append([node_to_index[tag] for tag in selected_tags])
        except KeyError as exc:
            raise ValueError(f"Gmsh element references unknown node tag {exc.args[0]}") from exc

    return vertices, np.asarray(faces, dtype=int)


def write_gmsh_mesh(
    mesh_path: str,
    vertices: np.ndarray,
    faces: np.ndarray,
    *,
    linearize: bool = False,
) -> None:
    """Write triangular mesh connectivity to an ASCII Gmsh 2.2 ``.msh`` file."""
    vertices = np.asarray(vertices, dtype=float)
    faces = np.asarray(faces, dtype=int)

    if vertices.ndim != 2 or vertices.shape[1] != 3:
        raise ValueError("vertices must have shape (n_vertices, 3)")
    if faces.ndim != 2 or faces.shape[1] < 3:
        raise ValueError("faces must have shape (n_faces, n_vertices_per_face) with at least 3 columns")

    element_rows = []
    for face in faces:
        active = face[face >= 0]
        if linearize:
            active = active[:3]
        element_type = _element_type_from_face(active)
        element_rows.append((element_type, active))

    path = Path(mesh_path)
    with path.open("w", encoding="utf-8") as stream:
        stream.write("$MeshFormat\n")
        stream.write("2.2 0 8\n")
        stream.write("$EndMeshFormat\n")
        stream.write("$Nodes\n")
        stream.write(f"{vertices.shape[0]}\n")
        for index, vertex in enumerate(vertices, start=1):
            stream.write(
                f"{index} "
                f"{vertex[0]:.17g} {vertex[1]:.17g} {vertex[2]:.17g}\n"
            )
        stream.write("$EndNodes\n")
        stream.write("$Elements\n")
        stream.write(f"{len(element_rows)}\n")
        for element_id, (element_type, active) in enumerate(element_rows, start=1):
            node_tags = " ".join(str(int(node) + 1) for node in active)
            stream.write(f"{element_id} {element_type} 0 {node_tags}\n")
        stream.write("$EndElements\n")


def _read_mesh_version(lines: Sequence[str]) -> str:
    start = _section_start(lines, "$MeshFormat")
    if start is None or start + 1 >= len(lines):
        raise ValueError("Gmsh file is missing a valid $MeshFormat section")
    return lines[start + 1].split()[0]


def _read_gmsh22(
    lines: Sequence[str],
    allowed_element_types: Iterable[int],
) -> Tuple[Dict[int, np.ndarray], list]:
    allowed = set(allowed_element_types)
    node_start = _require_section_start(lines, "$Nodes")
    n_nodes = int(lines[node_start + 1].strip())
    node_data = {}
    for line in lines[node_start + 2:node_start + 2 + n_nodes]:
        fields = line.split()
        node_data[int(fields[0])] = np.array([float(fields[1]), float(fields[2]), float(fields[3])])

    element_start = _require_section_start(lines, "$Elements")
    n_elements = int(lines[element_start + 1].strip())
    elements = []
    for line in lines[element_start + 2:element_start + 2 + n_elements]:
        fields = line.split()
        element_type = int(fields[1])
        if element_type not in allowed:
            continue
        n_tags = int(fields[2])
        node_tags = [int(value) for value in fields[3 + n_tags:]]
        elements.append((element_type, node_tags))

    return node_data, elements


def _read_gmsh4(
    lines: Sequence[str],
    allowed_element_types: Iterable[int],
) -> Tuple[Dict[int, np.ndarray], list]:
    allowed = set(allowed_element_types)

    node_tokens = _section_tokens(lines, "$Nodes", "$EndNodes")
    n_entity_blocks = int(next(node_tokens))
    int(next(node_tokens))
    int(next(node_tokens))
    int(next(node_tokens))

    node_data = {}
    for _ in range(n_entity_blocks):
        entity_dim = int(next(node_tokens))
        int(next(node_tokens))
        parametric = int(next(node_tokens))
        n_nodes_in_block = int(next(node_tokens))
        node_tags = [int(next(node_tokens)) for _ in range(n_nodes_in_block)]
        coordinate_size = 3 + (entity_dim if parametric else 0)
        for node_tag in node_tags:
            values = [float(next(node_tokens)) for _ in range(coordinate_size)]
            node_data[node_tag] = np.array(values[:3])

    elements = []
    element_lines = _section_lines(lines, "$Elements", "$EndElements")
    n_entity_blocks = int(element_lines[0].split()[0])
    cursor = 1
    for _ in range(n_entity_blocks):
        fields = element_lines[cursor].split()
        element_type = int(fields[2])
        n_elements_in_block = int(fields[3])
        cursor += 1
        for _ in range(n_elements_in_block):
            fields = element_lines[cursor].split()
            cursor += 1
            if element_type not in allowed:
                continue
            elements.append((element_type, [int(value) for value in fields[1:]]))

    return node_data, elements


def _section_start(lines: Sequence[str], name: str) -> Optional[int]:
    for index, line in enumerate(lines):
        if line.strip() == name:
            return index
    return None


def _require_section_start(lines: Sequence[str], name: str) -> int:
    start = _section_start(lines, name)
    if start is None:
        raise ValueError(f"Gmsh file is missing {name}")
    return start


def _section_tokens(lines: Sequence[str], start_name: str, end_name: str):
    start = _require_section_start(lines, start_name)
    tokens = []
    for line in lines[start + 1:]:
        if line.strip() == end_name:
            return iter(tokens)
        tokens.extend(line.split())
    raise ValueError(f"Gmsh file is missing {end_name}")


def _section_lines(lines: Sequence[str], start_name: str, end_name: str):
    start = _require_section_start(lines, start_name)
    section = []
    for line in lines[start + 1:]:
        if line.strip() == end_name:
            return section
        section.append(line)
    raise ValueError(f"Gmsh file is missing {end_name}")


def _element_type_from_face(face: np.ndarray) -> int:
    node_count = int(face.shape[0])
    element_type = _TRIANGLE_TYPE_BY_NODE_COUNT.get(node_count)
    if element_type is None:
        supported = ", ".join(str(count) for count in sorted(_TRIANGLE_TYPE_BY_NODE_COUNT))
        raise ValueError(
            f"Cannot write face with {node_count} nodes as a triangular Gmsh element. "
            f"Supported node counts are: {supported}"
        )
    return element_type
