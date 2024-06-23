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

def read_mesh_data(mesh_path: str) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """
    Read mesh data from a MAT file.

    Parameters
    ----------
    mesh_path : str
        The file path to the MAT file containing mesh data.

    Returns
    -------
    Tuple[Optional[np.ndarray], Optional[np.ndarray]]
        Vertices and faces data from the MAT file.

    Raises
    ------
    FileNotFoundError
        If the specified file does not exist.
    Exception
        If an error occurs during file reading.
    """
    if not os.path.exists(mesh_path):
        raise FileNotFoundError(f"File not found: {mesh_path}")

    try:
        mesh_mat = scipy.io.loadmat(mesh_path)
        key_list = list(mesh_mat.keys())
        vertices = mesh_mat[key_list[-1]]
        faces = mesh_mat[key_list[-2]] - 1  # Convert to zero-based indexing

        return vertices, faces
    except Exception as e:
        print(f"An error occurred while reading the mesh data: {e}")
        return None, None