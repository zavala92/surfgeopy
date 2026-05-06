"""High-level API for integrating functions over implicit surfaces."""

from dataclasses import dataclass, replace
from typing import Callable, Optional

import numpy as np

from .reference_quadrature import PULL_BACK_GAUSS
from .surf_integration import (
    DEFAULT_INTEGRATION_DEGREE,
    accumulate_surface_integrals,
    compute_surf_quadrature,
)
from .utils import read_mesh_data

__all__ = [
    "SurfaceMesh",
    "LevelSetSurface",
    "IntegrationConfig",
    "IntegrationResult",
    "DiagnosticIntegrationResult",
    "integrate",
    "integrate_with_diagnostics",
]


@dataclass(frozen=True)
class SurfaceMesh:
    """Triangulated reference mesh for an embedded surface.

    The mesh supplies the coarse reference triangles. The actual curved
    geometry is recovered by projecting interpolation nodes to a level set.
    """

    vertices: np.ndarray
    faces: np.ndarray

    def __post_init__(self) -> None:
        vertices = np.asarray(self.vertices, dtype=float)
        faces = np.asarray(self.faces, dtype=int)

        if vertices.ndim != 2 or vertices.shape[1] != 3:
            raise ValueError("vertices must have shape (n_vertices, 3)")
        if faces.ndim != 2 or faces.shape[1] < 3:
            raise ValueError("faces must have shape (n_faces, n_vertices_per_face) with at least 3 columns")

        object.__setattr__(self, "vertices", vertices)
        object.__setattr__(self, "faces", faces)

    @classmethod
    def from_mat(cls, mesh_path: str) -> "SurfaceMesh":
        """Load a MATLAB mesh file using surfgeopy's existing MAT convention."""
        vertices, faces = read_mesh_data(mesh_path)
        return cls(vertices, faces)

    @property
    def n_vertices(self) -> int:
        return self.vertices.shape[0]

    @property
    def n_faces(self) -> int:
        return self.faces.shape[0]


@dataclass(frozen=True)
class LevelSetSurface:
    """Implicit surface described by a reference mesh and level-set functions."""

    mesh: SurfaceMesh
    level_set: Callable[[np.ndarray], float]
    gradient: Callable[[np.ndarray], np.ndarray]


@dataclass(frozen=True)
class IntegrationConfig:
    """Numerical configuration for a surface integration run."""

    interpolation_degree: int
    lp_degree: float = float("inf")
    refinement_level: int = 0
    integration_degree: int = DEFAULT_INTEGRATION_DEGREE
    quadrature_rule: str = PULL_BACK_GAUSS

    def __post_init__(self) -> None:
        if self.interpolation_degree < 1:
            raise ValueError("interpolation_degree must be at least 1")
        if self.refinement_level < 0:
            raise ValueError("refinement_level must be non-negative")
        if self.integration_degree < 1:
            raise ValueError("integration_degree must be at least 1")


@dataclass(frozen=True)
class IntegrationResult:
    """Quadrature data and integrated values returned by :func:`integrate`."""

    values: np.ndarray
    points: np.ndarray
    weights: np.ndarray
    offsets: np.ndarray
    config: IntegrationConfig

    @property
    def total(self) -> float:
        """Return the total integral over the whole surface."""
        return float(np.sum(self.values))

    @property
    def n_quadrature_points(self) -> int:
        """Return the number of quadrature points used in this run."""
        return self.points.shape[0]

    def __float__(self) -> float:
        return self.total


@dataclass(frozen=True)
class DiagnosticIntegrationResult:
    """Error-estimation data returned by :func:`integrate_with_diagnostics`.

    The estimate is computed by comparing a base integration run with an
    enriched run using higher interpolation and/or quadrature degree.
    """

    base_result: IntegrationResult
    enriched_result: IntegrationResult
    absolute_error_estimate: float
    relative_error_estimate: float
    local_absolute_errors: np.ndarray
    recommended_config: IntegrationConfig
    absolute_tolerance: Optional[float] = None
    relative_tolerance: Optional[float] = 1.0e-8

    @property
    def total(self) -> float:
        """Return the enriched integral value."""
        return self.enriched_result.total

    @property
    def base_total(self) -> float:
        """Return the integral value from the base configuration."""
        return self.base_result.total

    @property
    def n_quadrature_points(self) -> int:
        """Return the number of quadrature points in the enriched run."""
        return self.enriched_result.n_quadrature_points

    @property
    def target_reached(self) -> Optional[bool]:
        """Return whether the requested accuracy target was reached.

        If no absolute or relative tolerance was supplied, ``None`` is returned.
        """
        checks = []
        if self.absolute_tolerance is not None:
            checks.append(self.absolute_error_estimate <= self.absolute_tolerance)
        if self.relative_tolerance is not None:
            checks.append(self.relative_error_estimate <= self.relative_tolerance)
        if not checks:
            return None
        return all(checks)

    @property
    def max_local_error(self) -> float:
        """Return the largest per-face difference between both runs."""
        if self.local_absolute_errors.size == 0:
            return 0.0
        return float(np.max(self.local_absolute_errors))

    def summary(self) -> str:
        """Return a compact text report for notebooks and logs."""
        status = self.target_reached
        if status is None:
            recommendation = "no tolerance requested"
        elif status:
            recommendation = "accuracy target reached"
        else:
            recommendation = "increase interpolation degree, quadrature degree, or refinement"

        return "\n".join(
            [
                f"Integral:              {self.total:.16g}",
                f"Base integral:         {self.base_total:.16g}",
                f"Estimated abs. error:  {self.absolute_error_estimate:.3e}",
                f"Estimated rel. error:  {self.relative_error_estimate:.3e}",
                f"Max local error:       {self.max_local_error:.3e}",
                f"Quadrature points:     {self.n_quadrature_points}",
                f"Interpolation degree:  {self.enriched_result.config.interpolation_degree}",
                f"Integration degree:    {self.enriched_result.config.integration_degree}",
                f"Recommendation:        {recommendation}",
            ]
        )


def integrate(
    surface: LevelSetSurface,
    integrand: Callable[[np.ndarray], float] = lambda _: 1.0,
    config: Optional[IntegrationConfig] = None,
) -> IntegrationResult:
    """Integrate a scalar function over an implicit surface.

    Parameters
    ----------
    surface
        Surface mesh plus implicit level-set representation.
    integrand
        Scalar function evaluated at physical quadrature points.
    config
        Numerical configuration. If omitted, a conservative default is used.

    Returns
    -------
    IntegrationResult
        Per-face values, quadrature points, weights, offsets, and total value.
    """
    if config is None:
        config = IntegrationConfig(interpolation_degree=6)

    points, weights, offsets = compute_surf_quadrature(
        surface.level_set,
        surface.gradient,
        surface.mesh.vertices,
        surface.mesh.faces,
        config.interpolation_degree,
        config.lp_degree,
        config.refinement_level,
        integrand,
        config.integration_degree,
        config.quadrature_rule,
    )
    values = accumulate_surface_integrals(points, weights, offsets, integrand)
    return IntegrationResult(values, points, weights, offsets, config)


def integrate_with_diagnostics(
    surface: LevelSetSurface,
    integrand: Callable[[np.ndarray], float] = lambda _: 1.0,
    config: Optional[IntegrationConfig] = None,
    *,
    interpolation_degree_step: int = 2,
    integration_degree_step: int = 2,
    absolute_tolerance: Optional[float] = None,
    relative_tolerance: Optional[float] = 1.0e-8,
) -> DiagnosticIntegrationResult:
    """Integrate a scalar function and estimate numerical accuracy.

    The function runs ``integrate`` twice: once with ``config`` and once with an
    enriched configuration. The absolute difference between the two totals is
    reported as an a posteriori error estimate.

    Parameters
    ----------
    surface
        Surface mesh plus implicit level-set representation.
    integrand
        Scalar function evaluated at physical quadrature points.
    config
        Base numerical configuration. If omitted, the same default as
        :func:`integrate` is used.
    interpolation_degree_step
        Increase applied to the interpolation degree for the enriched run.
    integration_degree_step
        Increase applied to the quadrature degree for the enriched run.
    absolute_tolerance
        Optional absolute error target.
    relative_tolerance
        Optional relative error target. The default is ``1e-8``.

    Returns
    -------
    DiagnosticIntegrationResult
        Base and enriched results, error estimates, local per-face differences,
        and a recommended next configuration.
    """
    if config is None:
        config = IntegrationConfig(interpolation_degree=6)
    if interpolation_degree_step < 0:
        raise ValueError("interpolation_degree_step must be non-negative")
    if integration_degree_step < 0:
        raise ValueError("integration_degree_step must be non-negative")
    if interpolation_degree_step == 0 and integration_degree_step == 0:
        raise ValueError("at least one enrichment step must be positive")

    base_result = integrate(surface, integrand, config)
    enriched_config = replace(
        config,
        interpolation_degree=config.interpolation_degree + interpolation_degree_step,
        integration_degree=config.integration_degree + integration_degree_step,
    )
    enriched_result = integrate(surface, integrand, enriched_config)

    absolute_error_estimate = abs(enriched_result.total - base_result.total)
    denominator = max(abs(enriched_result.total), np.finfo(float).tiny)
    relative_error_estimate = absolute_error_estimate / denominator

    if enriched_result.values.shape == base_result.values.shape:
        local_absolute_errors = np.abs(enriched_result.values - base_result.values)
    else:
        local_absolute_errors = np.array([], dtype=float)

    if (
        (absolute_tolerance is None or absolute_error_estimate <= absolute_tolerance)
        and (relative_tolerance is None or relative_error_estimate <= relative_tolerance)
    ):
        recommended_config = config
    else:
        recommended_config = replace(
            enriched_config,
            interpolation_degree=enriched_config.interpolation_degree + interpolation_degree_step,
            integration_degree=enriched_config.integration_degree + integration_degree_step,
        )

    return DiagnosticIntegrationResult(
        base_result=base_result,
        enriched_result=enriched_result,
        absolute_error_estimate=float(absolute_error_estimate),
        relative_error_estimate=float(relative_error_estimate),
        local_absolute_errors=local_absolute_errors,
        recommended_config=recommended_config,
        absolute_tolerance=absolute_tolerance,
        relative_tolerance=relative_tolerance,
    )
