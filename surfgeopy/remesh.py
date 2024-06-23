# remesh.py
"""
remesh.py
-------------
Deal with re-triangulation of existing meshes.
"""

import numpy as np
from .utils import *

__all__ = ['subdivide', 'faces_to_edges', 'unique_rows', 'hashable_rows', 'unique_ordered']

def subdivide(vertices, faces, face_index=None, vertex_attributes=None, return_index=False):
    """
    Subdivide a mesh into smaller triangles.
    Note that if `face_index` is passed, only those faces will be subdivided 
    and their neighbors won't be modified, making the mesh no longer "watertight."

    Parameters
    ----------
    vertices : (n, 3) float
        Vertices in space
    faces : (m, 3) int
        Indices of vertices which make up triangular faces
    face_index : array_like, optional
        Indices of faces to subdivide. If None, all faces of the mesh will be subdivided.
    vertex_attributes : dict, optional
        Contains (n, d) attribute data
    return_index : bool, optional
        If True, return index of original face for new faces

    Returns
    -------
    new_vertices : (q, 3) float
        Vertices in space
    new_faces : (p, 3) int
        Remeshed faces
    index_dict : dict, optional
        Only returned if `return_index` is True. 
        {index of original face : index of new faces}.
    """
    # Create a mask to select faces for subdivision
    face_mask = np.ones(len(faces), dtype=bool) if face_index is None else np.zeros(len(faces), dtype=bool)
    if face_index is not None:
        face_mask[face_index] = True

    # Select faces to be subdivided
    faces_subset = faces[face_mask]

    # Find unique edges of the selected faces
    edges = np.sort(faces_to_edges(faces_subset), axis=1)
    unique, inverse = unique_rows(edges)
    
    # Compute midpoints for unique edges
    mid = vertices[edges[unique]].mean(axis=1)
    mid_idx = inverse.reshape((-1, 3)) + len(vertices)

    # Create new faces using the original vertices and midpoints
    f = np.column_stack([
        faces_subset[:, 0], mid_idx[:, 0], mid_idx[:, 2],
        mid_idx[:, 0], faces_subset[:, 1], mid_idx[:, 1],
        mid_idx[:, 2], mid_idx[:, 1], faces_subset[:, 2],
        mid_idx[:, 0], mid_idx[:, 1], mid_idx[:, 2]
    ]).reshape((-1, 3))

    # Combine the old faces and new subdivided faces
    new_faces = np.vstack((faces[~face_mask], f))
    new_vertices = np.vstack((vertices, mid))

    # Handle vertex attributes if provided
    if vertex_attributes is not None:
        new_attributes = {}
        for key, values in vertex_attributes.items():
            attr_tris = values[faces_subset]
            attr_mid = np.vstack([attr_tris[:, g, :].mean(axis=1) for g in [[0, 1], [1, 2], [2, 0]]])
            attr_mid = attr_mid[unique]
            new_attributes[key] = np.vstack((values, attr_mid))
        return new_vertices, new_faces, new_attributes

    # Return index mapping of new faces to original faces if required
    if return_index:
        nonzero = np.nonzero(face_mask)[0]
        start = len(faces) - len(nonzero)
        stack = np.arange(start, start + len(f)).reshape((-1, 4))
        index_dict = {k: v for k, v in zip(nonzero, stack)}
        return new_vertices, new_faces, index_dict

    return new_vertices, new_faces

def faces_to_edges(faces, return_index=False):
    """
    Given a list of faces (n,3), return a list of edges (n*3,2)

    Parameters
    ----------
    faces : (n, 3) int
        Vertex indices representing faces

    Returns
    -------
    edges : (n*3, 2) int
        Vertex indices representing edges
    """
    faces = np.asanyarray(faces)

    # Each face has three edges
    edges = faces[:, [0, 1, 1, 2, 2, 0]].reshape((-1, 2))

    if return_index:
        # Create index array for edges
        face_index = np.tile(np.arange(len(faces)), (3, 1)).T.reshape(-1)
        return edges, face_index
    return edges

def unique_rows(data, digits=None, keep_order=False):
    """
    Returns indices of unique rows. It will return the first occurrence of a row that is duplicated:
    [[1,2], [3,4], [1,2]] will return [0,1]

    Parameters
    ----------
    data : (n, m) array
        Floating point data
    digits : int or None, optional
        Number of digits to consider for uniqueness

    Returns
    -------
    unique : (j,) int
        Index in data which is a unique row
    inverse : (n,) int
        Array to reconstruct original
        Example: data[unique][inverse] == data
    """
    # Convert rows to a hashable format
    rows = hashable_rows(data, digits=digits)
    if keep_order:
        return unique_ordered(rows, return_index=True, return_inverse=True)[1:]
    return np.unique(rows, return_index=True, return_inverse=True)[1:]

def hashable_rows(data, digits=None):
    """
    Convert array rows into a hashable format.

    Parameters
    ----------
    data : (n, m) array
        Input data
    digits : int or None, optional
        Number of digits to add to hash if data is floating point

    Returns
    -------
    hashable : (n,) array
        Custom data type which can be sorted or used as hash keys
    """
    if len(data) == 0:
        return np.array([])

    # Convert data to integer format based on precision
    as_int = float_to_int(data, digits=digits)
    if len(as_int.shape) == 1:
        return as_int

    # Use bitwise operations if the array is 2D and small enough
    if len(as_int.shape) == 2 and as_int.shape[1] <= 4:
        precision = int(np.floor(64 / as_int.shape[1]))
        if np.abs(as_int).max() < 2 ** (precision - 1):
            hashable = np.zeros(len(as_int), dtype=np.int64)
            for offset, column in enumerate(as_int.astype(np.int64).T):
                np.bitwise_xor(hashable, column << (offset * precision), out=hashable)
            return hashable

    # Use a custom data type for larger arrays
    dtype = np.dtype((np.void, as_int.dtype.itemsize * as_int.shape[1]))
    return np.ascontiguousarray(as_int).view(dtype).reshape(-1)

def unique_ordered(data, return_index=False, return_inverse=False):
    """
    Returns the same as np.unique, but ordered as per the first occurrence of the unique value in data.

    Parameters
    ----------
    data : array-like
        Input data
    return_index : bool, optional
        Return indices of unique values
    return_inverse : bool, optional
        Return the inverse of the unique array

    Returns
    -------
    unique : array-like
        The sorted unique values
    index : array-like, optional
        The indices of the unique values
    inverse : array-like, optional
        The indices to reconstruct the original data
    """
    unique, index, inverse = np.unique(data, return_index=True, return_inverse=True)
    order = index.argsort()

    result = [unique[order]]
    if return_index:
        result.append(index[order])
    if return_inverse:
        result.append(order.argsort()[inverse])

    return tuple(result) if len(result) > 1 else result[0]