from typing import Callable, Optional, Tuple

import numpy as np
from minterpy import MultiIndexSet, Grid, NewtonPolynomial
from minterpy.dds import dds

from .reference_quadrature import DEFAULT_QUADRATURE_RULE, make_reference_quadrature
from .remesh import subdivide
from .surface import ImplicitSurface, project_triangle_nodes, simplex_barycentric_coordinates
from .utils import (
    compute_norm, read_mesh_data, pushforward, _cross
)

__all__ = [
    'integration', 'accumulate_surface_integrals', 'compute_surf_quadrature',
    'compute_surf_geometry', 'quadrature_surf_tri', 'quadrature_split_surf_tri'
]

DEFAULT_INTEGRATION_DEGREE = 14
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
        quadrature_rule (Optional[str], optional): Reference quadrature rule,
            for example 'Pull_back_Gauss', 'Gauss_Legendre', or a 'ModePy_*'
            simplex rule. Defaults to None.

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
    quadrature_rule: str = DEFAULT_QUADRATURE_RULE
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
        quadrature_rule: str: Reference quadrature rule, for example
            'Pull_back_Gauss', 'Gauss_Legendre', or a 'ModePy_*' simplex rule.

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


def compute_surf_geometry(
    ls_function: Callable[[np.ndarray], float],
    ls_grad_func: Callable[[np.ndarray], np.ndarray],
    vertices: np.ndarray,
    faces: np.ndarray,
    interp_deg: int,
    lp_dgr: int,
    Refinement: int,
    deg_integration: int = 14,
    quadrature_rule: str = DEFAULT_QUADRATURE_RULE
) -> Tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    """Sample differential geometry quantities on the curved surface patches.

    The same Minterpy interpolant used by the integration routine is
    differentiated spectrally to obtain first and second parametric
    derivatives of the surface map.
    """
    if interp_deg < 2:
        raise ValueError("interp_deg must be at least 2 to compute curvature")

    surface = ImplicitSurface(ls_function, ls_grad_func)
    n_faces = faces.shape[0]
    nv_surf = faces.shape[1]
    offset = np.zeros(n_faces + 1, dtype=int)

    mi = MultiIndexSet.from_degree(spatial_dimension=2, poly_degree=interp_deg, lp_degree=lp_dgr)
    grid = Grid(mi)
    generating_points = pushforward(grid.unisolvent_nodes, duffy_transform=False)
    quad_ps = simplex_barycentric_coordinates(generating_points)
    reference_quadrature = make_reference_quadrature(deg_integration, quadrature_rule)

    points = []
    weights = []
    tangent_s = []
    tangent_t = []
    normal = []
    metric_tensor = []
    area_density = []
    second_fundamental_form = []
    mean_curvature = []
    gaussian_curvature = []

    for fun_id in range(n_faces):
        offset[fun_id] = len(points)

        n_elem = nv_surf - 1
        while faces[fun_id, n_elem] < 0:
            n_elem -= 1
        if n_elem < 2:
            continue

        for j in range(1, n_elem):
            lvids = [0, j, j + 1]
            local_vertices = vertices[faces[fun_id, lvids]]
            local_faces = np.array([[0, 1, 2]])
            for _ in range(Refinement):
                local_vertices, local_faces = subdivide(local_vertices, local_faces)

            for local_face in range(local_faces.shape[0]):
                pnts_p = project_triangle_nodes(
                    surface,
                    local_vertices[local_faces[local_face]],
                    quad_ps,
                )
                interpol_coeffs = np.squeeze(dds(pnts_p, grid.tree))
                newt_poly = NewtonPolynomial(mi, interpol_coeffs)
                ds_poly = newt_poly.diff([1, 0], backend="numba-par")
                dt_poly = newt_poly.diff([0, 1], backend="numba-par")
                dss_poly = newt_poly.diff([2, 0], backend="numba-par")
                dst_poly = newt_poly.diff([1, 1], backend="numba-par")
                dtt_poly = newt_poly.diff([0, 2], backend="numba-par")

                for qq in range(reference_quadrature.size):
                    evaluation_point = reference_quadrature.evaluation_points[qq]
                    evaluation_point = np.array([[evaluation_point[0], evaluation_point[1]]])
                    point = newt_poly(evaluation_point)[0]
                    p_s = ds_poly(evaluation_point)[0]
                    p_t = dt_poly(evaluation_point)[0]
                    p_ss = dss_poly(evaluation_point)[0]
                    p_st = dst_poly(evaluation_point)[0]
                    p_tt = dtt_poly(evaluation_point)[0]

                    cross = _cross(p_s, p_t)
                    jacobian = compute_norm(cross)
                    if jacobian <= np.finfo(float).eps:
                        unit_normal = np.full(3, np.nan)
                    else:
                        unit_normal = cross / jacobian

                    E = float(np.dot(p_s, p_s))
                    F = float(np.dot(p_s, p_t))
                    G = float(np.dot(p_t, p_t))
                    L = float(np.dot(p_ss, unit_normal))
                    M = float(np.dot(p_st, unit_normal))
                    N = float(np.dot(p_tt, unit_normal))
                    denominator = E * G - F * F
                    if denominator <= np.finfo(float).eps or not np.isfinite(denominator):
                        mean = np.nan
                        gaussian = np.nan
                    else:
                        mean = (E * N - 2.0 * F * M + G * L) / (2.0 * denominator)
                        gaussian = (L * N - M * M) / denominator

                    points.append(point)
                    weights.append(
                        reference_quadrature.weights[qq] *
                        jacobian *
                        reference_quadrature.weight_scale(qq)
                    )
                    tangent_s.append(p_s)
                    tangent_t.append(p_t)
                    normal.append(unit_normal)
                    metric_tensor.append([[E, F], [F, G]])
                    area_density.append(jacobian)
                    second_fundamental_form.append([[L, M], [M, N]])
                    mean_curvature.append(mean)
                    gaussian_curvature.append(gaussian)

    offset[n_faces] = len(points)
    return (
        np.asarray(points, dtype=float),
        np.asarray(weights, dtype=float),
        offset,
        np.asarray(tangent_s, dtype=float),
        np.asarray(tangent_t, dtype=float),
        np.asarray(normal, dtype=float),
        np.asarray(metric_tensor, dtype=float),
        np.asarray(area_density, dtype=float),
        np.asarray(second_fundamental_form, dtype=float),
        np.asarray(mean_curvature, dtype=float),
        np.asarray(gaussian_curvature, dtype=float),
    )

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
        quadrature_rule: str: Reference quadrature rule, for example
            'Pull_back_Gauss', 'Gauss_Legendre', or a 'ModePy_*' simplex rule.
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
        quadrature_rule: str: Reference quadrature rule, for example
            'Pull_back_Gauss', 'Gauss_Legendre', or a 'ModePy_*' simplex rule.
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
