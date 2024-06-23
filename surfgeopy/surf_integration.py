# Imports
import numpy as np
from minterpy import MultiIndexSet, Grid, NewtonPolynomial
from minterpy.dds import dds
from numba import njit
# Local imports
from .quadrature_points import quadrule_on_simplex
from .quadrature_points_gl import gauss_legendre_square
from .remesh import subdivide
from .utils import (
    SimpleImplicitSurfaceProjection, compute_norm, read_mesh_data,
    pushforward, pullback, _cross
)
from typing import Callable, Tuple, Optional

__all__ = ['integration', 'compute_surf_quadrature', 'quadrature_surf_tri', 'quadrature_split_surf_tri']

from typing import Callable, Optional
import numpy as np

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

    if deg_integration > 0:
        if quadrature_rule is None:
            quadrature_rule = 'Pull_back_Gauss'
        pnts, ws, offset = compute_surf_quadrature(
            ls_function, ls_grad_func, vertices, faces,
            interp_deg, lp_dgr, Refinement, fun_handle, deg_integration, quadrature_rule
        )
    else:
        pnts, ws, offset = compute_surf_quadrature(
            ls_function, ls_grad_func, vertices, faces,
            interp_deg, lp_dgr, Refinement, fun_handle
        )

    n_faces = faces.shape[0]
    fs = np.zeros(n_faces)

    for fun_id in range(n_faces):
        fs[fun_id] = np.sum(fun_handle(pnts[pid]) * ws[pid] for pid in range(offset[fun_id], offset[fun_id + 1]))

    return fs

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
    # Initialization
    index = 0
    n_faces = faces.shape[0]
    nv_surf = faces.shape[1]
    max_nv = max(1000000, n_faces * 6)
    pnts = np.zeros((max_nv, 3))
    ws = np.zeros(max_nv)
    offset = np.zeros(n_faces + 1, dtype=int)
    # Go through all the faces
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
        # Generate quadrature points
            if Refinement > 0:
                index = quadrature_split_surf_tri(
                    ls_function, ls_grad_func, pnts_tri, np.array([[0, 1, 2]]),
                    interp_deg, lp_dgr, Refinement, fun_handle, deg_integration,
                    quadrature_rule, pnts, ws, index
                )
            else:
                index = quadrature_surf_tri(
                    ls_function, ls_grad_func, pnts_tri, np.array([[0, 1, 2]]),
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
    n_faces = faces.shape[0]
    mi = MultiIndexSet.from_degree(spatial_dimension=2, poly_degree=interp_deg, lp_degree=lp_dgr)
    grid = Grid(mi)
    # Transform Chebyshev points from [-1,1]^2 to the reference simplex.
    generating_points = pushforward(grid.unisolvent_nodes, duffy_transform=False)
    quad_ps = np.array([[(1.0 - gp[0] - gp[1]), gp[0], gp[1]] for gp in generating_points])

    if quadrature_rule == 'Pull_back_Gauss':
        ws0, cs0 = quadrule_on_simplex(deg_integration)
        nqp = ws0.shape[0]
    # Transform quadrature points from the reference simplex to a unit square
        ksi = pullback(cs0, duffy_transform=False)
    else:
        ws0, cs0 = gauss_legendre_square(deg_integration)
        nqp = ws0.shape[0]
     # enlarge the size of quadrature points buffer if inadequate
    if index + n_faces * nqp > len(ws):
        n_new = 2 * len(ws) + n_faces * nqp
        ws.resize(n_new, refcheck=False)
        pnts.resize((n_new, 3), refcheck=False)
   
    for fun_id in range(n_faces):
        pnts_p = np.zeros((grid.unisolvent_nodes.shape[0], 3))
        for q, qp in enumerate(quad_ps):
            pnts_qq = (
                qp[0] * vertices[faces[fun_id, 0]] +
                qp[1] * vertices[faces[fun_id, 1]] +
                qp[2] * vertices[faces[fun_id, 2]]
            )
            pnts_p[q] = SimpleImplicitSurfaceProjection(ls_function, ls_grad_func, pnts_qq)

        interpol_coeffs = np.squeeze(dds(pnts_p, grid.tree))
        newt_poly = NewtonPolynomial(mi, interpol_coeffs)
        # compute partial derivatives with respect to "s"
        ds_poly = newt_poly.diff([1, 0], backend="numba-par")
         # compute partial derivatives with respect to "t"
        dt_poly = newt_poly.diff([0, 1], backend="numba-par")

        if quadrature_rule == 'Pull_back_Gauss':
            for qq in range(nqp):
                pnts[index] = newt_poly(np.array([ksi[qq, 0], ksi[qq, 1]]))
                # evaluate ∂_s at the quadrature points
                p_s = ds_poly(np.array([ksi[qq, 0], ksi[qq, 1]]))
                # evaluate ∂_t at the quadrature points
                p_t = dt_poly(np.array([ksi[qq, 0], ksi[qq, 1]]))
                # Compute ||∂_s x ∂_t||
                J = compute_norm(_cross(p_s, p_t))
                # Please use this in the case you are applying Duffy' transform
                #ws[index] = ws0[qq] * J * (4/(1-cs0[qq, 1]))
                ws[index] = ws0[qq] * J * (8 / np.sqrt((cs0[qq, 0] - cs0[qq, 1])**2 + 4 * (1 - cs0[qq, 0] - cs0[qq, 1])))
                index += 1
        else:
            for qq in range(nqp):
                pnts[index] = newt_poly(np.array([cs0[qq, 0], cs0[qq, 1]]))
                # evaluate ∂_s at the quadrature points
                p_s = ds_poly(np.array([cs0[qq, 0], cs0[qq, 1]]))
                # evaluate ∂_t at the quadrature points
                p_t = dt_poly(np.array([cs0[qq, 0], cs0[qq, 1]]))
                # Compute ||∂_s x ∂_t||
                J = compute_norm(_cross(p_s, p_t))
                ws[index] = ws0[qq] * J
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
