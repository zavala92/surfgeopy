# Imports
import numpy as np
import numba
from numba import njit
import scipy.io
import os
from typing import Callable, Optional, Tuple

__all__ = [
    'compute_norm', '_cross', 'max_edge_length', 'decimal_to_digits',
    'float_to_int', 'pushforward', 'pullback', 'SimpleImplicitSurfaceProjection', 'read_mesh_data'
]

_TYPE_MAP = [("f8", "i4"), ("f8", "i8")]
NB_OPTS = {"nogil": True}
# which works out to be 1e-13
TOL_ZERO = np.finfo(np.float64).resolution * 100
# how close to merge vertices
TOL_MERGE = 1e-8
_VERTEX_KEYS = ("vertices", "vertex", "verts", "xs", "points", "nodes", "coordinates")
_FACE_KEYS = ("faces", "face", "surfs", "triangles", "tris", "elements", "connectivity")

@njit(["float64(float64[:])"], **NB_OPTS)
def compute_norm(vec: np.ndarray) -> float:
    """
    Compute the Euclidean norm of a given vector.

    Parameters
    ----------
    vec : numpy.ndarray
        Input vector for which the Euclidean norm needs to be computed.

    Returns
    -------
    float
        Euclidean norm of the input vector.
    """
    sqnorm = np.float64(0.0)
    for i in range(len(vec)):
        sqnorm += vec[i] * vec[i]
    return np.sqrt(sqnorm)

@njit(["float64[:](float64[:], float64[:])"], **NB_OPTS)
def _cross(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """
    Compute the cross product of two vectors.

    Parameters
    ----------
    a, b : np.ndarray
        Input vectors.

    Returns
    -------
    np.ndarray
        Cross product of the input vectors.
    """
    return np.array([
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]
    ])

@njit(["float64(float64[:, :])"], **NB_OPTS)
def max_edge_length(xs: np.ndarray) -> float:
    """
    Compute the maximum edge length of a triangle defined by its vertices.

    Parameters
    ----------
    xs : numpy.ndarray
        Array of shape (3, N) representing N triangles, where each row contains the
        coordinates of the vertices of a triangle in 3D space.

    Returns
    -------
    float
        Maximum edge length among all the triangles.
    """
    return max(
        compute_norm(xs[0] - xs[1]),
        compute_norm(xs[1] - xs[2]),
        compute_norm(xs[2] - xs[0])
    )

def decimal_to_digits(decimal: float, min_digits: Optional[int] = None) -> int:
    """
    Return the number of digits to the first nonzero decimal.

    Parameters
    ----------
    decimal : float
        The decimal number to analyze.
    min_digits : Optional[int], default=None
        Minimum number of digits to return.

    Returns
    -------
    int
        Number of digits to the first nonzero decimal.
    """
    digits = abs(int(np.log10(decimal)))
    if min_digits is not None:
        digits = np.clip(digits, min_digits, 20)
    return digits

def float_to_int(data: np.ndarray, digits: Optional[int] = None, dtype=np.int32) -> np.ndarray:
    """
    Convert a numpy array of float/bool/int to integers.

    Parameters
    ----------
    data : np.ndarray
        Input data array.
    digits : Optional[int], default=None
        Precision for float conversion.
    dtype : np.dtype, default=np.int32
        Datatype for the result.

    Returns
    -------
    np.ndarray
        Data converted to integers.
    """
    data = np.asanyarray(data)
    if data.dtype.kind in 'ib' or data.size == 0:
        return data.astype(dtype)
    if data.dtype.kind != 'f':
        data = data.astype(np.float64)

    tol_merge = TOL_MERGE
    if digits is None:
        digits = decimal_to_digits(tol_merge)
    elif isinstance(digits, (float, np.float64)):
        digits = decimal_to_digits(digits)
    elif not isinstance(digits, (int, np.integer)):
        raise ValueError('Digits must be None, int, or float!')

    data_max = np.abs(data).max() * 10**digits
    dtype = np.int64 if data_max > 2**31 else np.int32
    as_int = np.round((data * 10 ** digits) - 1e-6).astype(dtype)
    return as_int

def pushforward(unisolvent_nodes: np.ndarray, duffy_transform: bool = False) -> np.ndarray:
    """
    Transform Chebyshev points from [-1, 1]^2 to a reference simplex.

    Parameters
    ----------
    unisolvent_nodes : np.ndarray
        Chebyshev points on the square.
    duffy_transform : bool, default=False
        Whether to apply Duffy's transform.

    Returns
    -------
    np.ndarray
        Transformed points on the simplex.
    """
    x, y = unisolvent_nodes[:, 0], unisolvent_nodes[:, 1]
    if duffy_transform:
        points_simplex_x = (1/4) * ((1 + x) * (1 - y))
        points_simplex_y = (1 + y) / 2
    else:
        points_simplex_x = (1 + x) * (3 - y) / 8
        points_simplex_y = (3 - x) * (y + 1) / 8

    return np.column_stack((points_simplex_x, points_simplex_y))

def pullback(qpoint_triangle: np.ndarray, duffy_transform: bool = False) -> np.ndarray:
    """
    Transform quadrature points from the reference simplex to a unit square.

    Parameters
    ----------
    qpoint_triangle : np.ndarray
        Quadrature points on the reference simplex.
    duffy_transform : bool, default=False
        Whether to apply Duffy's transform.

    Returns
    -------
    np.ndarray
        Transformed points on the [-1, 1]^2.
    """
    x, y = qpoint_triangle[:, 0], qpoint_triangle[:, 1]
    if duffy_transform:
        qpoint_square_x = (2 * x / (1 - y)) - 1
        qpoint_square_y = 2 * y - 1
    else:
        sqrt_term = np.sqrt((x - y) ** 2 + 4 * (1 - x - y))
        qpoint_square_x = 1 + (x - y) - sqrt_term
        qpoint_square_y = 1 - (x - y) - sqrt_term

    return np.column_stack((qpoint_square_x, qpoint_square_y))

def SimpleImplicitSurfaceProjection(
    phi: Callable[[np.ndarray], float],
    dphi: Callable[[np.ndarray], np.ndarray],
    x: np.ndarray,
    max_iter: int = 10
) -> np.ndarray:
    """
    Closest-point projection to surface given by an implicit function.

    Parameters
    ----------
    phi : Callable[[np.ndarray], float]
        Zero-levelset function.
    dphi : Callable[[np.ndarray], np.ndarray]
        Gradient of zero-levelset function.
    x : np.ndarray
        The point to be projected.
    max_iter : int, default=10
        Maximum number of iterations for the projection.

    Returns
    -------
    np.ndarray
        The projection point.
    """
    tol = 10 * np.finfo(np.float64).eps
    phi_v = phi(x)
    for _ in range(max_iter):
        grad_phi = dphi(x)
        grad_phi_norm = np.sum(grad_phi**2)
        normalize = phi_v / grad_phi_norm

        if np.sqrt(phi_v * normalize) < tol:
            break

        for j in range(len(x)):
            x[j] -= grad_phi[j] * normalize

        phi_v = phi(x)

    return x

def read_mesh_data(mesh_path: str) -> Tuple[np.ndarray, np.ndarray]:
    """
    Read mesh data from a MAT file.

    Parameters
    ----------
    mesh_path : str
        The file path to the MAT file containing mesh data.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        Vertices and faces data from the MAT file.

    Raises
    ------
    FileNotFoundError
        If the specified file does not exist.
    ValueError
        If the MAT file does not contain identifiable mesh arrays.
    """
    if not os.path.exists(mesh_path):
        raise FileNotFoundError(f"File not found: {mesh_path}")

    mesh_mat = scipy.io.loadmat(mesh_path)
    arrays = {
        key: value
        for key, value in mesh_mat.items()
        if not key.startswith("__") and isinstance(value, np.ndarray)
    }
    if not arrays:
        raise ValueError(f"No numeric array data found in MAT mesh file: {mesh_path}")

    vertex_key = _named_mesh_key(arrays, _VERTEX_KEYS)
    if vertex_key is None:
        vertex_key = _infer_vertex_key(arrays)
    vertices = _validate_vertices(arrays[vertex_key], vertex_key)

    face_key = _named_mesh_key(arrays, _FACE_KEYS)
    if face_key is None or face_key == vertex_key:
        face_key = _infer_face_key(arrays, vertices.shape[0], exclude={vertex_key})
    faces = _validate_faces(arrays[face_key], face_key, vertices.shape[0])

    return vertices, faces


def _named_mesh_key(arrays: dict, candidate_names: Tuple[str, ...]) -> Optional[str]:
    keys_by_lower_name = {key.lower(): key for key in arrays}
    for name in candidate_names:
        key = keys_by_lower_name.get(name)
        if key is not None:
            return key
    return None


def _infer_vertex_key(arrays: dict) -> str:
    candidates = [
        key for key, value in arrays.items()
        if _looks_like_vertices(value)
    ]
    if len(candidates) == 1:
        return candidates[0]
    if not candidates:
        raise ValueError(
            "Could not identify vertex coordinates in MAT file. "
            f"Expected one of {', '.join(_VERTEX_KEYS)} or a unique (n, 3) float array."
        )
    raise ValueError(
        "Could not uniquely identify vertex coordinates in MAT file. "
        f"Candidates: {', '.join(candidates)}"
    )


def _infer_face_key(arrays: dict, n_vertices: int, *, exclude: set) -> str:
    candidates = [
        key for key, value in arrays.items()
        if key not in exclude and _looks_like_faces(value, n_vertices)
    ]
    if len(candidates) == 1:
        return candidates[0]
    if not candidates:
        raise ValueError(
            "Could not identify face connectivity in MAT file. "
            f"Expected one of {', '.join(_FACE_KEYS)} or a unique integer connectivity array."
        )
    raise ValueError(
        "Could not uniquely identify face connectivity in MAT file. "
        f"Candidates: {', '.join(candidates)}"
    )


def _looks_like_vertices(value: np.ndarray) -> bool:
    array = _mesh_array(value)
    return (
        array.ndim == 2
        and array.shape[1] == 3
        and np.issubdtype(array.dtype, np.floating)
        and np.all(np.isfinite(array))
    )


def _looks_like_faces(value: np.ndarray, n_vertices: int) -> bool:
    array = _mesh_array(value)
    if array.ndim != 2 or array.shape[1] < 3 or not np.all(np.isfinite(array)):
        return False
    if not _is_integer_like(array):
        return False
    faces = np.rint(array).astype(int)
    if faces.size == 0:
        return True
    min_index = int(np.min(faces))
    max_index = int(np.max(faces))
    return (min_index >= 0 and max_index < n_vertices) or (min_index >= 1 and max_index <= n_vertices)


def _validate_vertices(value: np.ndarray, key: str) -> np.ndarray:
    vertices = _mesh_array(value).astype(float)
    if vertices.ndim != 2 or vertices.shape[1] != 3:
        raise ValueError(f"MAT variable {key!r} must have shape (n_vertices, 3)")
    if not np.all(np.isfinite(vertices)):
        raise ValueError(f"MAT variable {key!r} contains non-finite vertex coordinates")
    return vertices


def _validate_faces(value: np.ndarray, key: str, n_vertices: int) -> np.ndarray:
    faces_array = _mesh_array(value)
    if faces_array.ndim != 2 or faces_array.shape[1] < 3:
        raise ValueError(f"MAT variable {key!r} must have shape (n_faces, n_vertices_per_face)")
    if not np.all(np.isfinite(faces_array)):
        raise ValueError(f"MAT variable {key!r} contains non-finite face indices")
    if not _is_integer_like(faces_array):
        raise ValueError(f"MAT variable {key!r} must contain integer face indices")

    faces = np.rint(faces_array).astype(int)
    if faces.size == 0:
        return faces
    min_index = int(np.min(faces))
    max_index = int(np.max(faces))
    if min_index >= 1 and max_index <= n_vertices:
        faces = faces - 1
    elif min_index < 0 or max_index >= n_vertices:
        raise ValueError(
            f"MAT variable {key!r} references vertices outside the valid range "
            f"0..{n_vertices - 1} or 1..{n_vertices}"
        )
    return faces


def _is_integer_like(array: np.ndarray) -> bool:
    return np.all(np.isclose(array, np.rint(array), rtol=0.0, atol=1.0e-12))


def _mesh_array(value: np.ndarray) -> np.ndarray:
    array = np.asarray(value)
    if array.ndim > 2:
        array = np.squeeze(array)
    if array.ndim == 1:
        array = array.reshape(1, -1)
    return array
