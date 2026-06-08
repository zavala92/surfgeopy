"""Prototype singular and near-singular surface integration."""

from dataclasses import dataclass, replace
from typing import Callable, Optional, Tuple

import numpy as np
from minterpy import Grid, MultiIndexSet, NewtonPolynomial
from minterpy.dds import dds
from numpy.polynomial.chebyshev import chebvander2d
from numpy.polynomial.legendre import leggauss
from scipy import special
from scipy.optimize import minimize

from .api import LevelSetSurface
from .remesh import subdivide
from .surface import ImplicitSurface, project_triangle_nodes, simplex_barycentric_coordinates
from .utils import _cross, compute_norm, pushforward

__all__ = [
    "SingularIntegrationConfig",
    "SingularIntegralResult",
    "SingularDiagnosticResult",
    "ScreenedParametrixConfig",
    "ScreenedParametrixResult",
    "LaplaceSingleLayerOperator",
    "ScreenedLaplaceBeltramiParametrixOperator",
    "laplace_single_layer_potential",
    "screened_laplace_beltrami_parametrix",
]


@dataclass(frozen=True)
class SingularIntegrationConfig:
    """Numerical parameters for the prototype singular quadrature."""

    interpolation_degree: int = 8
    lp_degree: float = float("inf")
    refinement_level: int = 0
    regular_order: int = 18
    smooth_degree: int = 10
    moment_order: int = 48
    correction_order: int = 0
    near_threshold: float = 2.5
    singular_tolerance: float = 1.0e-10
    singular_model: str = "curvature"
    evaluation_strategy: str = "sspi"
    qbx_order: int = 12
    qbx_quadrature_order: int = 36
    qbx_radius_factor: float = 0.5
    qbx_radius_padding: float = 0.25
    qbx_diagnostic_order_step: int = 2

    def __post_init__(self) -> None:
        if self.interpolation_degree < 2:
            raise ValueError("interpolation_degree must be at least 2")
        if self.refinement_level < 0:
            raise ValueError("refinement_level must be non-negative")
        if self.regular_order < 1:
            raise ValueError("regular_order must be positive")
        if self.smooth_degree < 1:
            raise ValueError("smooth_degree must be positive")
        if self.moment_order < self.smooth_degree + 2:
            raise ValueError("moment_order should be at least smooth_degree + 2")
        if self.correction_order < 0:
            raise ValueError("correction_order must be non-negative")
        if self.near_threshold <= 0.0:
            raise ValueError("near_threshold must be positive")
        if self.singular_tolerance <= 0.0:
            raise ValueError("singular_tolerance must be positive")
        if self.singular_model not in {"metric", "curvature"}:
            raise ValueError("singular_model must be 'metric' or 'curvature'")
        if self.evaluation_strategy not in {"sspi", "hybrid_qbx"}:
            raise ValueError("evaluation_strategy must be 'sspi' or 'hybrid_qbx'")
        if self.qbx_order < 0:
            raise ValueError("qbx_order must be non-negative")
        if self.qbx_quadrature_order < 1:
            raise ValueError("qbx_quadrature_order must be positive")
        if self.qbx_radius_factor <= 0.0:
            raise ValueError("qbx_radius_factor must be positive")
        if self.qbx_radius_padding < 0.0:
            raise ValueError("qbx_radius_padding must be non-negative")
        if self.qbx_diagnostic_order_step < 0:
            raise ValueError("qbx_diagnostic_order_step must be non-negative")


@dataclass(frozen=True)
class SingularIntegralResult:
    """Potential values and panel classification counts."""

    values: np.ndarray
    near_panel_counts: np.ndarray
    singular_panel_counts: np.ndarray
    qbx_target_flags: Optional[np.ndarray] = None
    qbx_radii: Optional[np.ndarray] = None
    qbx_orders: Optional[np.ndarray] = None
    qbx_convergence_ratios: Optional[np.ndarray] = None
    qbx_error_estimates: Optional[np.ndarray] = None


@dataclass(frozen=True)
class SingularDiagnosticResult:
    """Diagnostic SSPI result from a base/enriched panel comparison."""

    values: np.ndarray
    base_values: np.ndarray
    enriched_values: np.ndarray
    absolute_error_estimates: np.ndarray
    tail_indicators: np.ndarray
    max_tail_indicators: np.ndarray
    near_panel_counts: np.ndarray
    singular_panel_counts: np.ndarray
    corrected_panel_counts: np.ndarray
    recommended_config: SingularIntegrationConfig
    target_tolerance: Optional[float] = 1.0e-8
    n_base_patches: int = 0
    n_enriched_patches_built: int = 0

    @property
    def enriched_patch_fraction(self) -> float:
        """Return the fraction of base patches that were enriched lazily."""
        if self.n_base_patches == 0:
            return 0.0
        return self.n_enriched_patches_built / self.n_base_patches


@dataclass(frozen=True)
class ScreenedParametrixConfig:
    """Numerical parameters for screened Laplace-Beltrami parametrix quadrature."""

    interpolation_degree: int = 8
    lp_degree: float = float("inf")
    refinement_level: int = 0
    regular_order: int = 18
    smooth_degree: int = 10
    moment_order: int = 48
    near_threshold: float = 2.5
    singular_tolerance: float = 1.0e-10
    singular_model: str = "curvature"
    curvature_corrected_distance: bool = True

    def __post_init__(self) -> None:
        if self.interpolation_degree < 2:
            raise ValueError("interpolation_degree must be at least 2")
        if self.refinement_level < 0:
            raise ValueError("refinement_level must be non-negative")
        if self.regular_order < 1:
            raise ValueError("regular_order must be positive")
        if self.smooth_degree < 1:
            raise ValueError("smooth_degree must be positive")
        if self.moment_order < self.smooth_degree + 2:
            raise ValueError("moment_order should be at least smooth_degree + 2")
        if self.near_threshold <= 0.0:
            raise ValueError("near_threshold must be positive")
        if self.singular_tolerance <= 0.0:
            raise ValueError("singular_tolerance must be positive")
        if self.singular_model not in {"metric", "curvature"}:
            raise ValueError("singular_model must be 'metric' or 'curvature'")


@dataclass(frozen=True)
class ScreenedParametrixResult:
    """Screened parametrix potential values and panel diagnostics."""

    values: np.ndarray
    near_panel_counts: np.ndarray
    singular_panel_counts: np.ndarray
    tail_indicators: np.ndarray


@dataclass(frozen=True)
class _PatchMap:
    polynomial: NewtonPolynomial
    ds_polynomial: NewtonPolynomial
    dt_polynomial: NewtonPolynomial
    dss_polynomial: NewtonPolynomial
    dst_polynomial: NewtonPolynomial
    dtt_polynomial: NewtonPolynomial
    diameter: float
    vertices: np.ndarray

    def evaluate(self, uv: np.ndarray) -> np.ndarray:
        return self.polynomial(_as_points(uv))

    def derivatives(self, uv: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        points = _as_points(uv)
        return self.ds_polynomial(points), self.dt_polynomial(points)

    def second_derivatives(self, uv: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        points = _as_points(uv)
        return (
            self.dss_polynomial(points),
            self.dst_polynomial(points),
            self.dtt_polynomial(points),
        )

    def jacobian(self, uv: np.ndarray) -> np.ndarray:
        ds, dt = self.derivatives(uv)
        return np.array([compute_norm(_cross(ds[i], dt[i])) for i in range(ds.shape[0])])

    def metric_at(self, uv: np.ndarray) -> np.ndarray:
        ds, dt = self.derivatives(np.asarray(uv, dtype=float).reshape(1, 2))
        s = ds[0]
        t = dt[0]
        return np.array([[np.dot(s, s), np.dot(s, t)], [np.dot(s, t), np.dot(t, t)]], dtype=float)

    def taylor_displacement(self, uv: np.ndarray, uv0: np.ndarray) -> np.ndarray:
        points = _as_points(uv)
        base = np.asarray(uv0, dtype=float).reshape(1, 2)
        diff = points - base
        ds, dt = self.derivatives(base)
        dss, dst, dtt = self.second_derivatives(base)

        s = ds[0]
        t = dt[0]
        ss = dss[0]
        st = dst[0]
        tt = dtt[0]
        xi = diff[:, 0]
        eta = diff[:, 1]
        linear = xi[:, None] * s + eta[:, None] * t
        quadratic = (
            0.5 * xi[:, None] * xi[:, None] * ss
            + xi[:, None] * eta[:, None] * st
            + 0.5 * eta[:, None] * eta[:, None] * tt
        )
        return linear + quadratic


@dataclass(frozen=True)
class _ClosestSurfaceData:
    patch: _PatchMap
    uv: np.ndarray
    point: np.ndarray
    normal: np.ndarray
    distance: float
    diameter: float
    signed_distance: float
    near_panel_count: int
    singular_panel_count: int
    use_qbx: bool


@dataclass(frozen=True)
class _QBXExpansionData:
    center: np.ndarray
    radius: float
    order: int
    diagnostic_order: int
    target_radius: float
    convergence_ratio: float


class LaplaceSingleLayerOperator:
    """Reusable prototype Laplace single-layer operator on a curved surface."""

    def __init__(
        self,
        surface: LevelSetSurface,
        config: Optional[SingularIntegrationConfig] = None,
    ) -> None:
        if config is None:
            config = SingularIntegrationConfig()
        self.surface = surface
        self.config = config
        self.patches = _build_patch_maps(surface, config)
        self._enriched_patch_cache = {}
        self._enriched_asset_cache = {}

    @property
    def n_patches(self) -> int:
        """Return the number of curved triangular patches in the operator."""
        return len(self.patches)

    def evaluate(
        self,
        targets: np.ndarray,
        density: Callable[[np.ndarray], float] = lambda _: 1.0,
    ) -> SingularIntegralResult:
        """Evaluate the potential at one or more target points."""
        if self.config.evaluation_strategy == "hybrid_qbx":
            return self.evaluate_hybrid_qbx(targets, density)

        target_points = _target_array(targets)
        values = np.zeros(target_points.shape[0], dtype=float)
        near_counts = np.zeros(target_points.shape[0], dtype=int)
        singular_counts = np.zeros(target_points.shape[0], dtype=int)

        for target_index, target in enumerate(target_points):
            total = 0.0
            for patch in self.patches:
                value, is_near, is_singular, _ = _evaluate_patch_contribution(
                    patch,
                    target,
                    density,
                    self.config,
                )
                near_counts[target_index] += int(is_near)
                singular_counts[target_index] += int(is_singular)
                total += value

            values[target_index] = total

        return SingularIntegralResult(values, near_counts, singular_counts)

    def evaluate_hybrid_qbx(
        self,
        targets: np.ndarray,
        density: Callable[[np.ndarray], float] = lambda _: 1.0,
    ) -> SingularIntegralResult:
        """Evaluate with SSPI away from the surface and QBX near the surface."""
        target_points = _target_array(targets)
        values = np.zeros(target_points.shape[0], dtype=float)
        near_counts = np.zeros(target_points.shape[0], dtype=int)
        singular_counts = np.zeros(target_points.shape[0], dtype=int)
        qbx_flags = np.zeros(target_points.shape[0], dtype=bool)
        qbx_radii = np.zeros(target_points.shape[0], dtype=float)
        qbx_orders = np.zeros(target_points.shape[0], dtype=int)
        qbx_ratios = np.zeros(target_points.shape[0], dtype=float)
        qbx_error_estimates = np.zeros(target_points.shape[0], dtype=float)

        for target_index, target in enumerate(target_points):
            closest = _closest_surface_data(self.patches, target, self.config)
            near_counts[target_index] = closest.near_panel_count
            singular_counts[target_index] = closest.singular_panel_count

            if closest.use_qbx:
                value, qbx_data, qbx_error_estimate = _qbx_single_layer_potential(
                    self.patches,
                    target,
                    density,
                    closest,
                    self.config,
                )
                values[target_index] = value
                qbx_flags[target_index] = True
                qbx_radii[target_index] = qbx_data.radius
                qbx_orders[target_index] = qbx_data.order
                qbx_ratios[target_index] = qbx_data.convergence_ratio
                qbx_error_estimates[target_index] = qbx_error_estimate
                continue

            total = 0.0
            for patch in self.patches:
                value, _, _, _ = _evaluate_patch_contribution(
                    patch,
                    target,
                    density,
                    self.config,
                )
                total += value
            values[target_index] = total

        return SingularIntegralResult(
            values,
            near_counts,
            singular_counts,
            qbx_flags,
            qbx_radii,
            qbx_orders,
            qbx_ratios,
            qbx_error_estimates,
        )

    def evaluate_with_diagnostics(
        self,
        targets: np.ndarray,
        density: Callable[[np.ndarray], float] = lambda _: 1.0,
        *,
        target_tolerance: Optional[float] = 1.0e-8,
        interpolation_degree_step: int = 2,
        regular_order_step: int = 2,
        smooth_degree_step: int = 2,
        moment_order_step: int = 16,
        correction_order_step: int = 8,
    ) -> SingularDiagnosticResult:
        """Evaluate with a base/enriched comparison on near and singular panels."""
        if target_tolerance is not None and target_tolerance <= 0.0:
            raise ValueError("target_tolerance must be positive or None")
        for name, step in (
            ("interpolation_degree_step", interpolation_degree_step),
            ("regular_order_step", regular_order_step),
            ("smooth_degree_step", smooth_degree_step),
            ("moment_order_step", moment_order_step),
            ("correction_order_step", correction_order_step),
        ):
            if step < 0:
                raise ValueError(f"{name} must be non-negative")

        enriched_config = replace(
            self.config,
            interpolation_degree=self.config.interpolation_degree + interpolation_degree_step,
            regular_order=self.config.regular_order + regular_order_step,
            smooth_degree=self.config.smooth_degree + smooth_degree_step,
            moment_order=self.config.moment_order + moment_order_step,
            correction_order=self.config.correction_order + correction_order_step,
        )
        if enriched_config.moment_order < enriched_config.smooth_degree + 2:
            enriched_config = replace(
                enriched_config,
                moment_order=enriched_config.smooth_degree + 2,
            )

        target_points = _target_array(targets)
        values = np.zeros(target_points.shape[0], dtype=float)
        base_values = np.zeros_like(values)
        enriched_values = np.zeros_like(values)
        absolute_error_estimates = np.zeros_like(values)
        tail_indicators = np.zeros_like(values)
        max_tail_indicators = np.zeros_like(values)
        near_counts = np.zeros(target_points.shape[0], dtype=int)
        singular_counts = np.zeros(target_points.shape[0], dtype=int)
        corrected_counts = np.zeros(target_points.shape[0], dtype=int)

        panel_tolerance = np.inf
        if target_tolerance is not None:
            panel_tolerance = target_tolerance / max(self.n_patches, 1)

        for target_index, target in enumerate(target_points):
            total = 0.0
            base_total = 0.0
            enriched_total = 0.0
            error_estimate = 0.0
            tail_indicator = 0.0
            max_tail_indicator = 0.0

            for patch_index, patch in enumerate(self.patches):
                base_value, is_near, is_singular, base_tail = _evaluate_patch_contribution(
                    patch,
                    target,
                    density,
                    self.config,
                )
                near_counts[target_index] += int(is_near)
                singular_counts[target_index] += int(is_singular)
                base_total += base_value
                tail_indicator += base_tail
                max_tail_indicator = max(max_tail_indicator, base_tail)

                if is_near or is_singular:
                    enriched_patch = self._enriched_patch(patch_index, enriched_config)
                    enriched_value, _, _, enriched_tail = _evaluate_patch_contribution(
                        enriched_patch,
                        target,
                        density,
                        enriched_config,
                    )
                    local_error = abs(enriched_value - base_value)
                    enriched_total += enriched_value
                    error_estimate += local_error

                    if local_error > panel_tolerance:
                        corrected_counts[target_index] += 1
                        total += enriched_value
                    else:
                        total += base_value
                else:
                    enriched_total += base_value
                    total += base_value

            values[target_index] = total
            base_values[target_index] = base_total
            enriched_values[target_index] = enriched_total
            absolute_error_estimates[target_index] = error_estimate
            tail_indicators[target_index] = tail_indicator
            max_tail_indicators[target_index] = max_tail_indicator

        recommended_config = self.config
        if np.any(corrected_counts > 0):
            recommended_config = enriched_config

        return SingularDiagnosticResult(
            values=values,
            base_values=base_values,
            enriched_values=enriched_values,
            absolute_error_estimates=absolute_error_estimates,
            tail_indicators=tail_indicators,
            max_tail_indicators=max_tail_indicators,
            near_panel_counts=near_counts,
            singular_panel_counts=singular_counts,
            corrected_panel_counts=corrected_counts,
            recommended_config=recommended_config,
            target_tolerance=target_tolerance,
            n_base_patches=self.n_patches,
            n_enriched_patches_built=len(self._enriched_patch_cache),
        )

    def _enriched_patch(
        self,
        patch_index: int,
        config: SingularIntegrationConfig,
    ) -> _PatchMap:
        key = (
            patch_index,
            config.interpolation_degree,
            config.lp_degree,
            config.refinement_level,
        )
        if key not in self._enriched_patch_cache:
            surface = ImplicitSurface(self.surface.level_set, self.surface.gradient)
            mi, grid, barycentric_points = self._patch_assets(config)
            self._enriched_patch_cache[key] = _build_patch_map(
                surface,
                self.patches[patch_index].vertices,
                config,
                mi=mi,
                grid=grid,
                barycentric_points=barycentric_points,
            )
        return self._enriched_patch_cache[key]

    def _patch_assets(
        self,
        config: SingularIntegrationConfig,
    ) -> Tuple[MultiIndexSet, Grid, np.ndarray]:
        key = (
            config.interpolation_degree,
            config.lp_degree,
        )
        if key not in self._enriched_asset_cache:
            self._enriched_asset_cache[key] = _patch_build_assets(config)
        return self._enriched_asset_cache[key]


class ScreenedLaplaceBeltramiParametrixOperator:
    """Prototype screened Laplace-Beltrami operator using a split Green kernel."""

    def __init__(
        self,
        surface: LevelSetSurface,
        alpha: float,
        config: Optional[ScreenedParametrixConfig] = None,
        smooth_remainder: Optional[Callable[[np.ndarray, np.ndarray], float]] = None,
    ) -> None:
        if alpha <= 0.0:
            raise ValueError("alpha must be positive")
        if config is None:
            config = ScreenedParametrixConfig()
        self.surface = surface
        self.alpha = float(alpha)
        self.config = config
        self.smooth_remainder = smooth_remainder
        self.patches = _build_patch_maps_from_screened_config(surface, config)

    @property
    def n_patches(self) -> int:
        """Return the number of curved triangular patches in the operator."""
        return len(self.patches)

    def evaluate(
        self,
        targets: np.ndarray,
        density: Callable[[np.ndarray], float] = lambda _: 1.0,
    ) -> ScreenedParametrixResult:
        """Evaluate the screened parametrix potential at one or more targets."""
        target_points = _target_array(targets)
        values = np.zeros(target_points.shape[0], dtype=float)
        near_counts = np.zeros(target_points.shape[0], dtype=int)
        singular_counts = np.zeros(target_points.shape[0], dtype=int)
        tail_indicators = np.zeros(target_points.shape[0], dtype=float)

        for target_index, target in enumerate(target_points):
            total = 0.0
            for patch in self.patches:
                value, is_near, is_singular, tail_indicator = (
                    _evaluate_screened_patch_contribution(
                        patch,
                        target,
                        density,
                        self.smooth_remainder,
                        self.alpha,
                        self.config,
                    )
                )
                near_counts[target_index] += int(is_near)
                singular_counts[target_index] += int(is_singular)
                tail_indicators[target_index] += tail_indicator
                total += value

            values[target_index] = total

        return ScreenedParametrixResult(
            values=values,
            near_panel_counts=near_counts,
            singular_panel_counts=singular_counts,
            tail_indicators=tail_indicators,
        )


def laplace_single_layer_potential(
    surface: LevelSetSurface,
    targets: np.ndarray,
    density: Callable[[np.ndarray], float] = lambda _: 1.0,
    config: Optional[SingularIntegrationConfig] = None,
) -> SingularIntegralResult:
    """Evaluate the Laplace single-layer potential with prototype product integration."""
    if config is None:
        config = SingularIntegrationConfig()
    return LaplaceSingleLayerOperator(surface, config).evaluate(targets, density)


def screened_laplace_beltrami_parametrix(
    surface: LevelSetSurface,
    targets: np.ndarray,
    alpha: float,
    density: Callable[[np.ndarray], float] = lambda _: 1.0,
    smooth_remainder: Optional[Callable[[np.ndarray, np.ndarray], float]] = None,
    config: Optional[ScreenedParametrixConfig] = None,
) -> ScreenedParametrixResult:
    """Evaluate the screened Laplace-Beltrami parametrix split prototype."""
    return ScreenedLaplaceBeltramiParametrixOperator(
        surface,
        alpha,
        config,
        smooth_remainder,
    ).evaluate(targets, density)


def _build_patch_maps(
    surface: LevelSetSurface,
    config: SingularIntegrationConfig,
) -> Tuple[_PatchMap, ...]:
    implicit_surface = ImplicitSurface(surface.level_set, surface.gradient)
    vertices = surface.mesh.vertices
    faces = surface.mesh.faces
    n_faces = faces.shape[0]
    nv_surf = faces.shape[1]

    mi, grid, barycentric_points = _patch_build_assets(config)
    patches = []

    for face_index in range(n_faces):
        n_elem = nv_surf - 1
        while faces[face_index, n_elem] < 0:
            n_elem -= 1
        if n_elem < 2:
            continue

        for local_index in range(1, n_elem):
            local_vertex_ids = [0, local_index, local_index + 1]
            local_vertices = vertices[faces[face_index, local_vertex_ids]]
            local_faces = np.array([[0, 1, 2]])
            for _ in range(config.refinement_level):
                local_vertices, local_faces = subdivide(local_vertices, local_faces)

            for local_face in local_faces:
                triangle_vertices = local_vertices[local_face]
                patches.append(
                    _build_patch_map(
                        implicit_surface,
                        triangle_vertices,
                        config,
                        mi=mi,
                        grid=grid,
                        barycentric_points=barycentric_points,
                    )
                )

    return tuple(patches)


def _build_patch_maps_from_screened_config(
    surface: LevelSetSurface,
    config: ScreenedParametrixConfig,
) -> Tuple[_PatchMap, ...]:
    singular_config = SingularIntegrationConfig(
        interpolation_degree=config.interpolation_degree,
        lp_degree=config.lp_degree,
        refinement_level=config.refinement_level,
        regular_order=config.regular_order,
        smooth_degree=config.smooth_degree,
        moment_order=config.moment_order,
        near_threshold=config.near_threshold,
        singular_tolerance=config.singular_tolerance,
        singular_model=config.singular_model,
    )
    return _build_patch_maps(surface, singular_config)


def _patch_build_assets(config: SingularIntegrationConfig) -> Tuple[MultiIndexSet, Grid, np.ndarray]:
    mi = MultiIndexSet.from_degree(
        spatial_dimension=2,
        poly_degree=config.interpolation_degree,
        lp_degree=config.lp_degree,
    )
    grid = Grid(mi)
    generating_points = pushforward(grid.unisolvent_nodes, duffy_transform=False)
    barycentric_points = simplex_barycentric_coordinates(generating_points)
    return mi, grid, barycentric_points


def _build_patch_map(
    surface: ImplicitSurface,
    triangle_vertices: np.ndarray,
    config: SingularIntegrationConfig,
    *,
    mi: Optional[MultiIndexSet] = None,
    grid: Optional[Grid] = None,
    barycentric_points: Optional[np.ndarray] = None,
) -> _PatchMap:
    if mi is None:
        mi = MultiIndexSet.from_degree(
            spatial_dimension=2,
            poly_degree=config.interpolation_degree,
            lp_degree=config.lp_degree,
        )
    if grid is None:
        grid = Grid(mi)
    if barycentric_points is None:
        generating_points = pushforward(grid.unisolvent_nodes, duffy_transform=False)
        barycentric_points = simplex_barycentric_coordinates(generating_points)

    projected_nodes = project_triangle_nodes(
        surface,
        triangle_vertices,
        barycentric_points,
    )
    coefficients = np.squeeze(dds(projected_nodes, grid.tree))
    polynomial = NewtonPolynomial(mi, coefficients)
    return _PatchMap(
        polynomial=polynomial,
        ds_polynomial=polynomial.diff([1, 0], backend="numba-par"),
        dt_polynomial=polynomial.diff([0, 1], backend="numba-par"),
        dss_polynomial=polynomial.diff([2, 0], backend="numba-par"),
        dst_polynomial=polynomial.diff([1, 1], backend="numba-par"),
        dtt_polynomial=polynomial.diff([0, 2], backend="numba-par"),
        diameter=_triangle_diameter(triangle_vertices),
        vertices=np.asarray(triangle_vertices, dtype=float),
    )


def _regular_integral(
    patch: _PatchMap,
    target: np.ndarray,
    density: Callable[[np.ndarray], float],
    order: int,
) -> float:
    nodes, weights = leggauss(order)
    uu, vv = np.meshgrid(nodes, nodes, indexing="ij")
    ww = np.outer(weights, weights)
    uv = np.column_stack((uu.ravel(), vv.ravel()))
    points = patch.evaluate(uv)
    jacobian = patch.jacobian(uv)

    total = 0.0
    for index in range(points.shape[0]):
        distance = compute_norm(points[index] - target)
        if distance <= np.finfo(float).tiny:
            continue
        total += (
            density(points[index])
            * jacobian[index]
            * ww.ravel()[index]
            / (4.0 * np.pi * distance)
        )
    return float(total)


def _evaluate_patch_contribution(
    patch: _PatchMap,
    target: np.ndarray,
    density: Callable[[np.ndarray], float],
    config: SingularIntegrationConfig,
) -> Tuple[float, bool, bool, float]:
    uv0, distance = _closest_square_parameter(patch, target)
    threshold = config.near_threshold * max(patch.diameter, np.finfo(float).eps)
    is_near = distance <= threshold
    is_singular = distance <= config.singular_tolerance * max(patch.diameter, 1.0)

    if is_near:
        value, tail_indicator = _product_integral_with_tail(
            patch,
            target,
            uv0,
            density,
            config,
        )
    else:
        value = _regular_integral(patch, target, density, config.regular_order)
        tail_indicator = 0.0

    return value, is_near, is_singular, tail_indicator


def _evaluate_screened_patch_contribution(
    patch: _PatchMap,
    target: np.ndarray,
    density: Callable[[np.ndarray], float],
    smooth_remainder: Optional[Callable[[np.ndarray, np.ndarray], float]],
    alpha: float,
    config: ScreenedParametrixConfig,
) -> Tuple[float, bool, bool, float]:
    uv0, distance = _closest_square_parameter(patch, target)
    threshold = config.near_threshold * max(patch.diameter, np.finfo(float).eps)
    is_near = distance <= threshold
    is_singular = distance <= config.singular_tolerance * max(patch.diameter, 1.0)

    if is_near:
        value, tail_indicator = _screened_product_integral_with_tail(
            patch,
            target,
            uv0,
            density,
            smooth_remainder,
            alpha,
            config,
        )
    else:
        value = _regular_screened_integral(
            patch,
            target,
            density,
            smooth_remainder,
            alpha,
            config,
        )
        tail_indicator = 0.0

    return value, is_near, is_singular, tail_indicator


def _regular_screened_integral(
    patch: _PatchMap,
    target: np.ndarray,
    density: Callable[[np.ndarray], float],
    smooth_remainder: Optional[Callable[[np.ndarray, np.ndarray], float]],
    alpha: float,
    config: ScreenedParametrixConfig,
) -> float:
    nodes, weights = leggauss(config.regular_order)
    uu, vv = np.meshgrid(nodes, nodes, indexing="ij")
    ww = np.outer(weights, weights).ravel()
    uv = np.column_stack((uu.ravel(), vv.ravel()))
    points = patch.evaluate(uv)
    jacobian = patch.jacobian(uv)
    distances = _screened_pair_distances(patch, uv, points, target, config)
    singular_kernel = _screened_singular_kernel(distances, alpha)

    total = 0.0
    for index, point in enumerate(points):
        kernel = singular_kernel[index]
        if smooth_remainder is not None:
            kernel += smooth_remainder(point, target)
        total += density(point) * jacobian[index] * kernel * ww[index]
    return float(total)


def _screened_product_integral_with_tail(
    patch: _PatchMap,
    target: np.ndarray,
    uv0: np.ndarray,
    density: Callable[[np.ndarray], float],
    smooth_remainder: Optional[Callable[[np.ndarray, np.ndarray], float]],
    alpha: float,
    config: ScreenedParametrixConfig,
) -> Tuple[float, float]:
    degree = config.smooth_degree
    cheb_nodes = np.cos(np.pi * np.arange(degree + 1) / degree)
    uu, vv = np.meshgrid(cheb_nodes, cheb_nodes, indexing="ij")
    uv = np.column_stack((uu.ravel(), vv.ravel()))

    points = patch.evaluate(uv)
    jacobian = patch.jacobian(uv)
    smooth_values = np.array(
        [
            density(points[index]) * jacobian[index]
            for index in range(points.shape[0])
        ],
        dtype=float,
    )

    vandermonde = chebvander2d(uv[:, 0], uv[:, 1], [degree, degree])
    coefficients = np.linalg.solve(vandermonde, smooth_values)
    moments = _screened_chebyshev_moments(
        patch,
        uv0,
        target,
        config.singular_model,
        config.curvature_corrected_distance,
        degree,
        config.moment_order,
        alpha,
    )
    singular_value = float(np.dot(coefficients, moments))
    tail_indicator = _chebyshev_tail_indicator(coefficients, moments, degree)

    smooth_value = 0.0
    if smooth_remainder is not None:
        smooth_value = _smooth_remainder_integral(
            patch,
            target,
            density,
            smooth_remainder,
            config.regular_order,
        )

    return singular_value + smooth_value, tail_indicator


def _smooth_remainder_integral(
    patch: _PatchMap,
    target: np.ndarray,
    density: Callable[[np.ndarray], float],
    smooth_remainder: Callable[[np.ndarray, np.ndarray], float],
    order: int,
) -> float:
    nodes, weights = leggauss(order)
    uu, vv = np.meshgrid(nodes, nodes, indexing="ij")
    ww = np.outer(weights, weights).ravel()
    uv = np.column_stack((uu.ravel(), vv.ravel()))
    points = patch.evaluate(uv)
    jacobian = patch.jacobian(uv)

    total = 0.0
    for index, point in enumerate(points):
        total += (
            density(point)
            * jacobian[index]
            * smooth_remainder(point, target)
            * ww[index]
        )
    return float(total)


def _product_integral(
    patch: _PatchMap,
    target: np.ndarray,
    uv0: np.ndarray,
    density: Callable[[np.ndarray], float],
    config: SingularIntegrationConfig,
) -> float:
    value, _ = _product_integral_with_tail(patch, target, uv0, density, config)
    return value


def _product_integral_with_tail(
    patch: _PatchMap,
    target: np.ndarray,
    uv0: np.ndarray,
    density: Callable[[np.ndarray], float],
    config: SingularIntegrationConfig,
) -> Tuple[float, float]:
    degree = config.smooth_degree
    cheb_nodes = np.cos(np.pi * np.arange(degree + 1) / degree)
    uu, vv = np.meshgrid(cheb_nodes, cheb_nodes, indexing="ij")
    uv = np.column_stack((uu.ravel(), vv.ravel()))

    points = patch.evaluate(uv)
    jacobian = patch.jacobian(uv)
    metric = patch.metric_at(uv0)
    residual = patch.evaluate(uv0)[0] - target
    model_distance = _model_distance(patch, uv, uv0, metric, residual, config.singular_model)
    exact_distance = np.linalg.norm(points - target, axis=1)
    ratio = np.ones_like(exact_distance)
    mask = exact_distance > 100.0 * np.finfo(float).eps
    ratio[mask] = model_distance[mask] / exact_distance[mask]
    smooth_values = np.array(
        [
            density(points[index]) * jacobian[index] * ratio[index] / (4.0 * np.pi)
            for index in range(points.shape[0])
        ]
    )

    vandermonde = chebvander2d(uv[:, 0], uv[:, 1], [degree, degree])
    coefficients = np.linalg.solve(vandermonde, smooth_values)
    moments = _singular_chebyshev_moments(
        patch,
        uv0,
        metric,
        residual,
        config.singular_model,
        degree,
        config.moment_order,
    )
    product_value = float(np.dot(coefficients, moments))
    tail_indicator = _chebyshev_tail_indicator(coefficients, moments, degree)

    if config.correction_order <= 0:
        return product_value, tail_indicator

    correction = _product_residual_correction(
        patch,
        target,
        uv0,
        metric,
        residual,
        config.singular_model,
        density,
        degree,
        coefficients,
        config.correction_order,
    )
    return product_value + correction, tail_indicator


def _product_residual_correction(
    patch: _PatchMap,
    target: np.ndarray,
    uv0: np.ndarray,
    metric: np.ndarray,
    residual: np.ndarray,
    singular_model: str,
    density: Callable[[np.ndarray], float],
    degree: int,
    coefficients: np.ndarray,
    order: int,
) -> float:
    nodes, weights = leggauss(order)
    uu, vv = np.meshgrid(nodes, nodes, indexing="ij")
    ww = np.outer(weights, weights).ravel()
    uv = np.column_stack((uu.ravel(), vv.ravel()))

    points = patch.evaluate(uv)
    jacobian = patch.jacobian(uv)
    exact_distance = np.linalg.norm(points - target, axis=1)
    model_distance = _model_distance(patch, uv, uv0, metric, residual, singular_model)
    exact_distance = np.maximum(exact_distance, np.finfo(float).tiny)
    model_distance = np.maximum(model_distance, np.finfo(float).tiny)

    basis = chebvander2d(uv[:, 0], uv[:, 1], [degree, degree])
    smooth_approximation = basis @ coefficients
    true_integrand = np.array(
        [
            density(points[index]) * jacobian[index] / (4.0 * np.pi * exact_distance[index])
            for index in range(points.shape[0])
        ]
    )
    model_integrand = smooth_approximation / model_distance
    return float(np.dot(ww, true_integrand - model_integrand))


def _closest_surface_data(
    patches: Tuple[_PatchMap, ...],
    target: np.ndarray,
    config: SingularIntegrationConfig,
) -> _ClosestSurfaceData:
    best_patch = patches[0]
    best_uv = np.zeros(2, dtype=float)
    best_distance = np.inf
    near_count = 0
    singular_count = 0

    for patch in patches:
        uv, distance = _closest_square_parameter(patch, target)
        threshold = config.near_threshold * max(patch.diameter, np.finfo(float).eps)
        near_count += int(distance <= threshold)
        singular_count += int(
            distance <= config.singular_tolerance * max(patch.diameter, 1.0)
        )
        if distance < best_distance:
            best_patch = patch
            best_uv = uv
            best_distance = distance

    closest_point = best_patch.evaluate(best_uv)[0]
    normal = _patch_normal(best_patch, best_uv)
    signed_distance = float(np.dot(target - closest_point, normal))
    qbx_threshold = config.near_threshold * max(best_patch.diameter, np.finfo(float).eps)
    return _ClosestSurfaceData(
        patch=best_patch,
        uv=best_uv,
        point=closest_point,
        normal=normal,
        distance=float(best_distance),
        diameter=best_patch.diameter,
        signed_distance=signed_distance,
        near_panel_count=near_count,
        singular_panel_count=singular_count,
        use_qbx=best_distance <= qbx_threshold,
    )


def _qbx_single_layer_potential(
    patches: Tuple[_PatchMap, ...],
    target: np.ndarray,
    density: Callable[[np.ndarray], float],
    closest: _ClosestSurfaceData,
    config: SingularIntegrationConfig,
) -> Tuple[float, _QBXExpansionData, float]:
    qbx_data = _select_qbx_expansion(target, closest, config)
    total = _qbx_single_layer_potential_at_order(
        patches,
        qbx_data.center,
        target,
        density,
        qbx_data.order,
        config.qbx_quadrature_order,
    )

    error_estimate = 0.0
    if qbx_data.diagnostic_order < qbx_data.order:
        diagnostic_total = _qbx_single_layer_potential_at_order(
            patches,
            qbx_data.center,
            target,
            density,
            qbx_data.diagnostic_order,
            config.qbx_quadrature_order,
        )
        error_estimate = abs(total - diagnostic_total)

    return float(total), qbx_data, float(error_estimate)


def _select_qbx_expansion(
    target: np.ndarray,
    closest: _ClosestSurfaceData,
    config: SingularIntegrationConfig,
) -> _QBXExpansionData:
    diameter = max(closest.diameter, np.finfo(float).eps)
    target_gap = abs(closest.signed_distance)
    radius = max(
        config.qbx_radius_factor * diameter,
        target_gap + config.qbx_radius_padding * diameter,
        100.0 * np.finfo(float).eps,
    )
    side = 1.0 if closest.signed_distance >= 0.0 else -1.0
    center = closest.point + side * radius * closest.normal
    target_radius = compute_norm(target - center)

    if target_radius > radius:
        radius = target_radius + max(
            config.qbx_radius_padding * diameter,
            100.0 * np.finfo(float).eps,
        )
        center = closest.point + side * radius * closest.normal
        target_radius = compute_norm(target - center)

    order = int(config.qbx_order)
    diagnostic_order = max(0, order - int(config.qbx_diagnostic_order_step))
    convergence_ratio = target_radius / max(radius, np.finfo(float).eps)
    return _QBXExpansionData(
        center=center,
        radius=float(radius),
        order=order,
        diagnostic_order=diagnostic_order,
        target_radius=float(target_radius),
        convergence_ratio=float(convergence_ratio),
    )


def _qbx_single_layer_potential_at_order(
    patches: Tuple[_PatchMap, ...],
    center: np.ndarray,
    target: np.ndarray,
    density: Callable[[np.ndarray], float],
    expansion_order: int,
    quadrature_order: int,
) -> float:
    total = 0.0
    for patch in patches:
        total += _qbx_patch_contribution(
            patch,
            center,
            target,
            density,
            expansion_order,
            quadrature_order,
        )
    return float(total)


def _qbx_patch_contribution(
    patch: _PatchMap,
    center: np.ndarray,
    target: np.ndarray,
    density: Callable[[np.ndarray], float],
    expansion_order: int,
    quadrature_order: int,
) -> float:
    nodes, weights = leggauss(quadrature_order)
    uu, vv = np.meshgrid(nodes, nodes, indexing="ij")
    ww = np.outer(weights, weights).ravel()
    uv = np.column_stack((uu.ravel(), vv.ravel()))
    points = patch.evaluate(uv)
    jacobian = patch.jacobian(uv)
    kernel = _qbx_laplace_kernel(center, target, points, expansion_order)
    values = np.array([density(point) for point in points], dtype=float)
    return float(np.dot(ww, values * jacobian * kernel))


def _qbx_laplace_kernel(
    center: np.ndarray,
    target: np.ndarray,
    sources: np.ndarray,
    expansion_order: int,
) -> np.ndarray:
    target_vector = np.asarray(target, dtype=float) - np.asarray(center, dtype=float)
    target_radius = compute_norm(target_vector)
    source_vectors = np.asarray(sources, dtype=float) - np.asarray(center, dtype=float).reshape(
        1,
        3,
    )
    source_radii = np.linalg.norm(source_vectors, axis=1)
    source_radii = np.maximum(source_radii, np.finfo(float).tiny)

    series = 1.0 / source_radii
    if expansion_order == 0 or target_radius <= np.finfo(float).eps:
        return series / (4.0 * np.pi)

    cos_angle = (source_vectors @ target_vector) / (source_radii * target_radius)
    cos_angle = np.clip(cos_angle, -1.0, 1.0)
    p_nm2 = np.ones_like(cos_angle)
    p_nm1 = cos_angle
    radius_power = target_radius
    series = series + radius_power * p_nm1 / (source_radii ** 2)

    for degree in range(2, expansion_order + 1):
        p_n = ((2 * degree - 1) * cos_angle * p_nm1 - (degree - 1) * p_nm2) / degree
        radius_power *= target_radius
        series = series + radius_power * p_n / (source_radii ** (degree + 1))
        p_nm2 = p_nm1
        p_nm1 = p_n

    return series / (4.0 * np.pi)


def _patch_normal(patch: _PatchMap, uv: np.ndarray) -> np.ndarray:
    ds, dt = patch.derivatives(np.asarray(uv, dtype=float).reshape(1, 2))
    normal = _cross(ds[0], dt[0])
    normal_norm = compute_norm(normal)
    if normal_norm <= np.finfo(float).tiny:
        return np.array([0.0, 0.0, 1.0])
    return normal / normal_norm


def _chebyshev_tail_indicator(
    coefficients: np.ndarray,
    moments: np.ndarray,
    degree: int,
) -> float:
    tail_width = max(1, min(degree, degree // 3))
    tail_start = degree - tail_width + 1
    mode_i, mode_j = np.indices((degree + 1, degree + 1))
    tail_mask = ((mode_i >= tail_start) | (mode_j >= tail_start)).ravel()
    weighted_tail = np.asarray(coefficients)[tail_mask] * np.asarray(moments)[tail_mask]
    return float(np.sum(np.abs(weighted_tail)))


def _singular_chebyshev_moments(
    patch: _PatchMap,
    uv0: np.ndarray,
    metric: np.ndarray,
    residual: np.ndarray,
    singular_model: str,
    degree: int,
    order: int,
) -> np.ndarray:
    if np.all(np.abs(uv0) <= 1.0 + 1.0e-13):
        return _duffy_chebyshev_moments(
            patch,
            uv0,
            metric,
            residual,
            singular_model,
            degree,
            order,
        )

    nodes, weights = leggauss(order)
    uu, vv = np.meshgrid(nodes, nodes, indexing="ij")
    ww = np.outer(weights, weights).ravel()
    uv = np.column_stack((uu.ravel(), vv.ravel()))
    denominator = _model_distance(patch, uv, uv0, metric, residual, singular_model)
    denominator = np.maximum(denominator, np.finfo(float).tiny)
    basis = chebvander2d(uv[:, 0], uv[:, 1], [degree, degree])
    return basis.T @ (ww / denominator)


def _duffy_chebyshev_moments(
    patch: _PatchMap,
    uv0: np.ndarray,
    metric: np.ndarray,
    residual: np.ndarray,
    singular_model: str,
    degree: int,
    order: int,
) -> np.ndarray:
    nodes, weights = leggauss(order)
    r = 0.5 * (nodes + 1.0)
    wr = 0.5 * weights
    s = r
    ws = wr

    moments = np.zeros((degree + 1) * (degree + 1), dtype=float)
    x_lengths = [uv0[0] + 1.0, 1.0 - uv0[0]]
    y_lengths = [uv0[1] + 1.0, 1.0 - uv0[1]]
    x_signs = [-1.0, 1.0]
    y_signs = [-1.0, 1.0]

    for x_length, x_sign in zip(x_lengths, x_signs):
        if x_length <= 1.0e-14:
            continue
        for y_length, y_sign in zip(y_lengths, y_signs):
            if y_length <= 1.0e-14:
                continue
            moments += _duffy_rectangle_moments(
                patch,
                uv0,
                metric,
                residual,
                singular_model,
                degree,
                r,
                wr,
                s,
                ws,
                x_length,
                y_length,
                x_sign,
                y_sign,
            )

    return moments


def _duffy_rectangle_moments(
    patch: _PatchMap,
    uv0: np.ndarray,
    metric: np.ndarray,
    residual: np.ndarray,
    singular_model: str,
    degree: int,
    r_nodes: np.ndarray,
    r_weights: np.ndarray,
    s_nodes: np.ndarray,
    s_weights: np.ndarray,
    x_length: float,
    y_length: float,
    x_sign: float,
    y_sign: float,
) -> np.ndarray:
    rr, ss = np.meshgrid(r_nodes, s_nodes, indexing="ij")
    ww = np.outer(r_weights, s_weights).ravel()
    rr = rr.ravel()
    ss = ss.ravel()

    uv_first = np.column_stack(
        (
            uv0[0] + x_sign * x_length * rr,
            uv0[1] + y_sign * y_length * rr * ss,
        )
    )
    jacobian_first = x_length * y_length * rr

    uv_second = np.column_stack(
        (
            uv0[0] + x_sign * x_length * rr * ss,
            uv0[1] + y_sign * y_length * rr,
        )
    )
    jacobian_second = x_length * y_length * rr

    moments = np.zeros((degree + 1) * (degree + 1), dtype=float)
    for uv, jacobian in ((uv_first, jacobian_first), (uv_second, jacobian_second)):
        denominator = _model_distance(patch, uv, uv0, metric, residual, singular_model)
        denominator = np.maximum(denominator, np.finfo(float).tiny)
        basis = chebvander2d(uv[:, 0], uv[:, 1], [degree, degree])
        moments += basis.T @ (ww * jacobian / denominator)

    return moments


def _screened_chebyshev_moments(
    patch: _PatchMap,
    uv0: np.ndarray,
    target: np.ndarray,
    singular_model: str,
    curvature_corrected_distance: bool,
    degree: int,
    order: int,
    alpha: float,
) -> np.ndarray:
    if np.all(np.abs(uv0) <= 1.0 + 1.0e-13):
        return _duffy_screened_chebyshev_moments(
            patch,
            uv0,
            target,
            singular_model,
            curvature_corrected_distance,
            degree,
            order,
            alpha,
        )

    nodes, weights = leggauss(order)
    uu, vv = np.meshgrid(nodes, nodes, indexing="ij")
    ww = np.outer(weights, weights).ravel()
    uv = np.column_stack((uu.ravel(), vv.ravel()))
    distances = _screened_model_distance(
        patch,
        uv,
        uv0,
        target,
        singular_model,
        curvature_corrected_distance,
    )
    basis = chebvander2d(uv[:, 0], uv[:, 1], [degree, degree])
    return basis.T @ (ww * _screened_singular_kernel(distances, alpha))


def _duffy_screened_chebyshev_moments(
    patch: _PatchMap,
    uv0: np.ndarray,
    target: np.ndarray,
    singular_model: str,
    curvature_corrected_distance: bool,
    degree: int,
    order: int,
    alpha: float,
) -> np.ndarray:
    nodes, weights = leggauss(order)
    r = 0.5 * (nodes + 1.0)
    wr = 0.5 * weights
    s = r
    ws = wr

    moments = np.zeros((degree + 1) * (degree + 1), dtype=float)
    x_lengths = [uv0[0] + 1.0, 1.0 - uv0[0]]
    y_lengths = [uv0[1] + 1.0, 1.0 - uv0[1]]
    x_signs = [-1.0, 1.0]
    y_signs = [-1.0, 1.0]

    for x_length, x_sign in zip(x_lengths, x_signs):
        if x_length <= 1.0e-14:
            continue
        for y_length, y_sign in zip(y_lengths, y_signs):
            if y_length <= 1.0e-14:
                continue
            moments += _duffy_rectangle_screened_moments(
                patch,
                uv0,
                target,
                singular_model,
                curvature_corrected_distance,
                degree,
                r,
                wr,
                s,
                ws,
                x_length,
                y_length,
                x_sign,
                y_sign,
                alpha,
            )

    return moments


def _duffy_rectangle_screened_moments(
    patch: _PatchMap,
    uv0: np.ndarray,
    target: np.ndarray,
    singular_model: str,
    curvature_corrected_distance: bool,
    degree: int,
    r_nodes: np.ndarray,
    r_weights: np.ndarray,
    s_nodes: np.ndarray,
    s_weights: np.ndarray,
    x_length: float,
    y_length: float,
    x_sign: float,
    y_sign: float,
    alpha: float,
) -> np.ndarray:
    rr, ss = np.meshgrid(r_nodes, s_nodes, indexing="ij")
    ww = np.outer(r_weights, s_weights).ravel()
    rr = rr.ravel()
    ss = ss.ravel()

    uv_first = np.column_stack(
        (
            uv0[0] + x_sign * x_length * rr,
            uv0[1] + y_sign * y_length * rr * ss,
        )
    )
    jacobian_first = x_length * y_length * rr

    uv_second = np.column_stack(
        (
            uv0[0] + x_sign * x_length * rr * ss,
            uv0[1] + y_sign * y_length * rr,
        )
    )
    jacobian_second = x_length * y_length * rr

    moments = np.zeros((degree + 1) * (degree + 1), dtype=float)
    for uv, jacobian in ((uv_first, jacobian_first), (uv_second, jacobian_second)):
        distances = _screened_model_distance(
            patch,
            uv,
            uv0,
            target,
            singular_model,
            curvature_corrected_distance,
        )
        basis = chebvander2d(uv[:, 0], uv[:, 1], [degree, degree])
        kernel = _screened_singular_kernel(distances, alpha)
        moments += basis.T @ (ww * jacobian * kernel)

    return moments


def _closest_square_parameter(patch: _PatchMap, target: np.ndarray) -> Tuple[np.ndarray, float]:
    coarse = np.linspace(-1.0, 1.0, 5)
    uu, vv = np.meshgrid(coarse, coarse, indexing="ij")
    candidates = np.column_stack((uu.ravel(), vv.ravel()))
    points = patch.evaluate(candidates)
    distances = np.sum((points - target) ** 2, axis=1)
    start = candidates[int(np.argmin(distances))]

    def objective(uv: np.ndarray) -> float:
        point = patch.evaluate(np.asarray(uv, dtype=float).reshape(1, 2))[0]
        residual = point - target
        return float(np.dot(residual, residual))

    result = minimize(objective, start, method="L-BFGS-B", bounds=[(-1.0, 1.0), (-1.0, 1.0)])
    uv0 = np.asarray(result.x if result.success else start, dtype=float)
    distance = np.sqrt(objective(uv0))
    return uv0, float(distance)


def _metric_distance(
    uv: np.ndarray,
    uv0: np.ndarray,
    metric: np.ndarray,
    delta: float,
) -> np.ndarray:
    diff = np.asarray(uv, dtype=float) - np.asarray(uv0, dtype=float).reshape(1, 2)
    quadratic = np.einsum("ni,ij,nj->n", diff, metric, diff)
    return np.sqrt(np.maximum(quadratic + delta * delta, 0.0))


def _model_distance(
    patch: _PatchMap,
    uv: np.ndarray,
    uv0: np.ndarray,
    metric: np.ndarray,
    residual: np.ndarray,
    singular_model: str,
) -> np.ndarray:
    if singular_model == "metric":
        return _metric_distance(uv, uv0, metric, compute_norm(residual))
    if singular_model == "curvature":
        displacement = patch.taylor_displacement(uv, uv0)
        local_vectors = displacement + np.asarray(residual, dtype=float).reshape(1, 3)
        quartic = np.sum(local_vectors * local_vectors, axis=1)
        return np.sqrt(np.maximum(quartic, 0.0))
    raise ValueError("singular_model must be 'metric' or 'curvature'")


def _screened_model_distance(
    patch: _PatchMap,
    uv: np.ndarray,
    uv0: np.ndarray,
    target: np.ndarray,
    singular_model: str,
    curvature_corrected_distance: bool,
) -> np.ndarray:
    points = patch.evaluate(uv)
    closest_point = patch.evaluate(uv0)[0]
    metric = patch.metric_at(uv0)
    residual = closest_point - target
    chord_distance = _model_distance(patch, uv, uv0, metric, residual, singular_model)
    if not curvature_corrected_distance:
        return np.maximum(chord_distance, np.finfo(float).tiny)

    gaussian_curvature = _patch_gaussian_curvature_at(patch, uv0)
    return _curvature_corrected_distance(chord_distance, gaussian_curvature)


def _screened_pair_distances(
    patch: _PatchMap,
    uv: np.ndarray,
    points: np.ndarray,
    target: np.ndarray,
    config: ScreenedParametrixConfig,
) -> np.ndarray:
    chord_distance = np.linalg.norm(points - target.reshape(1, 3), axis=1)
    chord_distance = np.maximum(chord_distance, np.finfo(float).tiny)
    if not config.curvature_corrected_distance:
        return chord_distance

    gaussian_curvature = _patch_gaussian_curvature_values(patch, uv)
    return _curvature_corrected_distance(chord_distance, gaussian_curvature)


def _screened_singular_kernel(distance: np.ndarray, alpha: float) -> np.ndarray:
    distance = np.maximum(np.asarray(distance, dtype=float), np.finfo(float).tiny)
    return special.k0(np.sqrt(alpha) * distance) / (2.0 * np.pi)


def _curvature_corrected_distance(
    chord_distance: np.ndarray,
    gaussian_curvature,
) -> np.ndarray:
    chord = np.asarray(chord_distance, dtype=float)
    curvature = np.asarray(gaussian_curvature, dtype=float)
    corrected = chord * (
        1.0
        + curvature * chord * chord / 24.0
        + 3.0 * curvature * curvature * chord**4 / 640.0
    )
    return np.maximum(corrected, np.finfo(float).tiny)


def _patch_gaussian_curvature_at(patch: _PatchMap, uv: np.ndarray) -> float:
    return float(_patch_gaussian_curvature_values(patch, np.asarray(uv).reshape(1, 2))[0])


def _patch_gaussian_curvature_values(patch: _PatchMap, uv: np.ndarray) -> np.ndarray:
    points = _as_points(uv)
    ds, dt = patch.derivatives(points)
    dss, dst, dtt = patch.second_derivatives(points)
    normal = np.cross(ds, dt)
    normal_norm = np.linalg.norm(normal, axis=1)
    safe_norm = np.maximum(normal_norm, np.finfo(float).tiny)
    normal = normal / safe_norm[:, None]

    first_eg = np.sum(ds * ds, axis=1) * np.sum(dt * dt, axis=1)
    first_eg -= np.sum(ds * dt, axis=1) ** 2
    second_eg = np.sum(normal * dss, axis=1) * np.sum(normal * dtt, axis=1)
    second_eg -= np.sum(normal * dst, axis=1) ** 2
    denominator = np.maximum(np.abs(first_eg), np.finfo(float).tiny)
    return second_eg / denominator


def _triangle_diameter(vertices: np.ndarray) -> float:
    return float(
        max(
            compute_norm(vertices[0] - vertices[1]),
            compute_norm(vertices[1] - vertices[2]),
            compute_norm(vertices[2] - vertices[0]),
        )
    )


def _as_points(uv: np.ndarray) -> np.ndarray:
    points = np.asarray(uv, dtype=float)
    if points.ndim == 1:
        return points.reshape(1, 2)
    return points


def _target_array(targets: np.ndarray) -> np.ndarray:
    target_points = np.asarray(targets, dtype=float)
    if target_points.ndim == 1:
        target_points = target_points.reshape(1, -1)
    if target_points.ndim != 2 or target_points.shape[1] != 3:
        raise ValueError("targets must have shape (n_targets, 3)")
    return target_points
