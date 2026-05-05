from typing import Callable, Optional, Tuple

import numpy as np
from minterpy import MultiIndexSet, Grid, NewtonPolynomial
from minterpy.dds import dds

from .reference_quadrature import PULL_BACK_GAUSS, make_reference_quadrature
from .remesh import subdivide
from .surface import ImplicitSurface, project_triangle_nodes, simplex_barycentric_coordinates
from .utils import (
    compute_norm, read_mesh_data, pushforward, _cross
)

__all__ = [
    'integration', 'accumulate_surface_integrals', 'compute_surf_quadrature',
    'quadrature_surf_tri', 'quadrature_split_surf_tri'
]

DEFAULT_INTEGRATION_DEGREE = 14
DEFAULT_QUADRATURE_RULE = PULL_BACK_GAUSS

def integration(
    ls_function: Callable[[np.ndarray], float],
    ls_grad_func: Callable[[np.ndarray], np.ndarray],
    mesh: str,
    interp_deg: int,
    lp_dgr: int,
    Refinement: int,
    fun_handle: Callable[[np.ndarray], float] = lambda _: 1.0,
    deg_integration: int = -1,
    quadrature_rule: Optional[str] = None
) -> np.ndarray:
    """
    Compute the integration of a function over curved triangles.

    Args:
        ls_function (Callable[[np.ndarray], float]): Zero-levelset function.
        ls_grad_func (Callable[[np.ndarray], np.ndarray]): Gradient of the zero-levelset function.
        mesh (str): The file path to the MAT file containing mesh data.
        interp_deg (int): Interpolation degree.
        lp_dgr (int): The l_p-norm used to define the polynomial degree.
        Refinement (int): Refinement level.
        fun_handle (Callable[[np.ndarray], float], optional): Function to be evaluated on each quadrature point. Defaults to a constant function.
        deg_integration (int, optional): Degree of integration. Defaults to -1 (use default configuration).
        quadrature_rule (Optional[str], optional): Quadrature rule type. Can be 'Gauss_Legendre' or 'Gauss_Simplex'. Defaults to None.

    Returns:
        np.ndarray: Integration values for each curved triangle.
    """
    vertices, faces = read_mesh_data(mesh)

    if deg_integration <= 0:
        deg_integration = DEFAULT_INTEGRATION_DEGREE
    if quadrature_rule is None:
        quadrature_rule = DEFAULT_QUADRATURE_RULE

    pnts, ws, offset = compute_surf_quadrature(
        ls_function, ls_grad_func, vertices, faces,
        interp_deg, lp_dgr, Refinement, fun_handle, deg_integration, quadrature_rule
    )
    fs = accumulate_surface_integrals(pnts, ws, offset, fun_handle)

    return fs


def accumulate_surface_integrals(
    pnts: np.ndarray,
    ws: np.ndarray,
    offset: np.ndarray,
    fun_handle: Callable[[np.ndarray], float],
) -> np.ndarray:
    """Accumulate pointwise quadrature data into one integral per face."""
    values = np.zeros(len(offset) - 1)
    for fun_id in range(len(values)):
        values[fun_id] = sum(
            fun_handle(pnts[pid]) * ws[pid]
            for pid in range(offset[fun_id], offset[fun_id + 1])
        )
    return values

def compute_surf_quadrature(
    ls_function: Callable[[np.ndarray], float],
    ls_grad_func: Callable[[np.ndarray], np.ndarray],
    vertices: np.ndarray,
    faces: np.ndarray,
    interp_deg: int,
    lp_dgr: int,
    Refinement: int,
    fun_handle: Callable[[np.ndarray], float],
    deg_integration: int = 14,
    quadrature_rule: str = 'Pull_back_Gauss'
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute quadrature points and weights on curved triangles.

    Args:
        ls_function: Callable[[np.ndarray], float]: Zero-levelset function.
        ls_grad_func: Callable[[np.ndarray], np.ndarray]: Gradient of zero-levelset function.
        vertices: np.ndarray: Array of vertex coordinates.
        faces: np.ndarray: Array of face connectivity.
        interp_deg: int: Interpolation degree.
        lp_dgr: int: :math:`l_p`-norm, which is used to define the polynomial degree.
        Refinement: int: Refinement level.
        fun_handle: Callable[[np.ndarray], float]: Function to be evaluated on each quadrature point.
        deg_integration: int: Degree of integration (default: 14).
        quadrature_rule: str: Quadrature rule type ('Gauss_Legendre' or 'Gauss_Simplex').

    Returns:
        Tuple[np.ndarray, np.ndarray, np.ndarray]: Quadrature points, weights, and offset array.
    """
    surface = ImplicitSurface(ls_function, ls_grad_func)

    index = 0
    n_faces = faces.shape[0]
    nv_surf = faces.shape[1]
    max_nv = max(1000000, n_faces * 6)
    pnts = np.zeros((max_nv, 3))
    ws = np.zeros(max_nv)
    offset = np.zeros(n_faces + 1, dtype=int)

    for fun_id in range(n_faces):
        offset[fun_id] = index

        n_elem = nv_surf - 1
        while faces[fun_id, n_elem] < 0:
            n_elem -= 1
        if n_elem < 2:
            continue
        # Split each element into several curved triangles
        for j in range(1, n_elem):
            lvids = [0, j, j + 1]
            pnts_tri = vertices[faces[fun_id, lvids]]

            if Refinement > 0:
                index = quadrature_split_surf_tri(
                    surface.level_set, surface.gradient, pnts_tri, np.array([[0, 1, 2]]),
                    interp_deg, lp_dgr, Refinement, fun_handle, deg_integration,
                    quadrature_rule, pnts, ws, index
                )
            else:
                index = quadrature_surf_tri(
                    surface.level_set, surface.gradient, pnts_tri, np.array([[0, 1, 2]]),
                    interp_deg, lp_dgr, fun_handle, deg_integration,
                    quadrature_rule, pnts, ws, index
                )

    pnts = pnts[:index]
    ws = ws[:index]
    offset[n_faces] = index
    return pnts, ws, offset

def quadrature_surf_tri(
    ls_function: Callable[[np.ndarray], float],
    ls_grad_func: Callable[[np.ndarray], np.ndarray],
    vertices: np.ndarray,
    faces: np.ndarray,
    interp_deg: int,
    lp_dgr: int,
    fun_handle: Callable[[np.ndarray], float],
    deg_integration: int,
    quadrature_rule: str,
    pnts: np.ndarray,
    ws: np.ndarray,
    index: int
) -> int:
    """
    For a mixed mesh, find the cell integration of the test function f.

    Args:
        ls_function: Callable[[np.ndarray], float]: Zero-levelset function.
        ls_grad_func: Callable[[np.ndarray], np.ndarray]: Gradient of zero-levelset function.
        vertices: np.ndarray: Array of vertex coordinates.
        faces: np.ndarray: Array of face connectivity.
        interp_deg: int: Interpolation degree.
        lp_dgr: int: :math:`l_p`-norm, which is used to define the polynomial degree.
        fun_handle: Callable[[np.ndarray], float]: Function to be evaluated on each quadrature point.
        deg_integration: int: Degree of integration.
        quadrature_rule: str: Quadrature rule type ('Gauss_Legendre' or 'Gauss_Simplex').
        pnts: np.ndarray: Quadrature points array.
        ws: np.ndarray: Quadrature weights array.
        index: int: Current index in the arrays.

    Returns:
        int: Updated index value.
    """
    surface = ImplicitSurface(ls_function, ls_grad_func)
    n_faces = faces.shape[0]
    mi = MultiIndexSet.from_degree(spatial_dimension=2, poly_degree=interp_deg, lp_degree=lp_dgr)
    grid = Grid(mi)
    generating_points = pushforward(grid.unisolvent_nodes, duffy_transform=False)
    quad_ps = simplex_barycentric_coordinates(generating_points)
    reference_quadrature = make_reference_quadrature(deg_integration, quadrature_rule)
    nqp = reference_quadrature.size

    if index + n_faces * nqp > len(ws):
        n_new = 2 * len(ws) + n_faces * nqp
        ws.resize(n_new, refcheck=False)
        pnts.resize((n_new, 3), refcheck=False)
   
    for fun_id in range(n_faces):
        pnts_p = project_triangle_nodes(surface, vertices[faces[fun_id]], quad_ps)

        interpol_coeffs = np.squeeze(dds(pnts_p, grid.tree))
        newt_poly = NewtonPolynomial(mi, interpol_coeffs)
        ds_poly = newt_poly.diff([1, 0], backend="numba-par")
        dt_poly = newt_poly.diff([0, 1], backend="numba-par")

        for qq in range(reference_quadrature.size):
            evaluation_point = reference_quadrature.evaluation_points[qq]
            evaluation_point = np.array([[evaluation_point[0], evaluation_point[1]]])
            pnts[index] = newt_poly(evaluation_point)[0]
            p_s = ds_poly(evaluation_point)[0]
            p_t = dt_poly(evaluation_point)[0]
            J = compute_norm(_cross(p_s, p_t))
            ws[index] = (
                reference_quadrature.weights[qq] *
                J *
                reference_quadrature.weight_scale(qq)
            )
            index += 1

    return index

def quadrature_split_surf_tri(
    ls_function: Callable[[np.ndarray], float],
    ls_grad_func: Callable[[np.ndarray], np.ndarray],
    vertices: np.ndarray,
    faces: np.ndarray,
    interp_deg: int,
    lp_dgr: int,
    Refinement: int,
    fun_handle: Callable[[np.ndarray], float],
    deg_integration: int,
    quadrature_rule: str,
    pnts: np.ndarray,
    ws: np.ndarray,
    index: int
) -> int:
    """
    For a mixed mesh, find the cell integration of the test function f.

    Args:
        ls_function: Callable[[np.ndarray], float]: Zero-levelset function.
        ls_grad_func: Callable[[np.ndarray], np.ndarray]: Gradient of the zero-levelset function.
        vertices: np.ndarray: Array of vertex coordinates.
        faces: np.ndarray: Array of face connectivity.
        interp_deg: int: Interpolation degree.
        lp_dgr: int: :math:`l_p`-norm, which is used to define the polynomial degree.
        deg_integration: int: Degree of integration.
        quadrature_rule: str: Quadrature rule type ('Gauss_Legendre' or 'Gauss_Simplex').
        pnts: np.ndarray: Quadrature points array.
        ws: np.ndarray: Quadrature weights array.
        index: int: Current index in the arrays.

    Returns:
        int: Updated index value.
    """
    for _ in range(Refinement):
        vertices, faces = subdivide(vertices, faces)

    index = quadrature_surf_tri(
        ls_function, ls_grad_func, vertices, faces, interp_deg,
        lp_dgr, fun_handle, deg_integration, quadrature_rule, pnts, ws, index
    )

    return index
