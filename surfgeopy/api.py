"""High-level API for integrating functions over implicit surfaces."""

from dataclasses import dataclass, replace
from typing import Callable, Optional, Sequence, Tuple

import numpy as np

from .mesh_io import read_gmsh_mesh, write_gmsh_mesh
from .reference_quadrature import DEFAULT_QUADRATURE_RULE
from .remesh import subdivide, subdivide_conforming
from .surf_integration import (
    DEFAULT_INTEGRATION_DEGREE,
    compute_surf_geometry,
    compute_surf_quadrature,
)
from .utils import read_mesh_data

__all__ = [
    "SurfaceMesh",
    "LevelSetSurface",
    "IntegrationConfig",
    "IntegrationResult",
    "CompiledIntegrationScheme",
    "DiagnosticIntegrationResult",
    "IndicatorRefinementIteration",
    "IndicatorRefinementResult",
    "PatchGeometry",
    "SurfaceGeometryResult",
    "compile_integration",
    "integrate",
    "integrate_with_diagnostics",
    "refine_by_indicator",
    "surface_geometry",
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

    @classmethod
    def icosphere(
        cls,
        *,
        refinement_level: int = 1,
        radius: float = 1.0,
    ) -> "SurfaceMesh":
        """Create a triangular icosphere reference mesh.

        This is a convenient built-in mesh for examples, tests, and first
        experiments on spherical level sets. Each refinement step splits every
        triangle and projects the new vertices back to the sphere.
        """
        if refinement_level < 0:
            raise ValueError("refinement_level must be non-negative")
        if radius <= 0.0:
            raise ValueError("radius must be positive")

        golden_ratio = (1.0 + np.sqrt(5.0)) / 2.0
        vertices = np.array(
            [
                [-1.0, golden_ratio, 0.0],
                [1.0, golden_ratio, 0.0],
                [-1.0, -golden_ratio, 0.0],
                [1.0, -golden_ratio, 0.0],
                [0.0, -1.0, golden_ratio],
                [0.0, 1.0, golden_ratio],
                [0.0, -1.0, -golden_ratio],
                [0.0, 1.0, -golden_ratio],
                [golden_ratio, 0.0, -1.0],
                [golden_ratio, 0.0, 1.0],
                [-golden_ratio, 0.0, -1.0],
                [-golden_ratio, 0.0, 1.0],
            ],
            dtype=float,
        )
        faces = np.array(
            [
                [0, 11, 5],
                [0, 5, 1],
                [0, 1, 7],
                [0, 7, 10],
                [0, 10, 11],
                [1, 5, 9],
                [5, 11, 4],
                [11, 10, 2],
                [10, 7, 6],
                [7, 1, 8],
                [3, 9, 4],
                [3, 4, 2],
                [3, 2, 6],
                [3, 6, 8],
                [3, 8, 9],
                [4, 9, 5],
                [2, 4, 11],
                [6, 2, 10],
                [8, 6, 7],
                [9, 8, 1],
            ],
            dtype=int,
        )

        vertices = _normalize_to_radius(vertices, radius)
        for _ in range(refinement_level):
            vertices, faces = subdivide(vertices, faces)
            vertices = _normalize_to_radius(vertices, radius)

        return cls(vertices, faces)

    @classmethod
    def from_gmsh(
        cls,
        mesh_path: str,
        *,
        preserve_order: bool = False,
        element_types: Optional[Sequence[int]] = None,
    ) -> "SurfaceMesh":
        """Load triangular surface elements from an ASCII Gmsh ``.msh`` file.

        High-order triangles are linearized to their corner nodes by default.
        Use ``preserve_order=True`` to keep all high-order element nodes in
        ``faces`` for inspection, conversion, or external processing.
        """
        vertices, faces = read_gmsh_mesh(
            mesh_path,
            preserve_order=preserve_order,
            element_types=element_types,
        )
        return cls(vertices, faces)

    def to_gmsh(self, mesh_path: str, *, linearize: bool = False) -> None:
        """Write this triangular mesh to an ASCII Gmsh 2.2 ``.msh`` file."""
        write_gmsh_mesh(mesh_path, self.vertices, self.faces, linearize=linearize)

    @property
    def n_vertices(self) -> int:
        return self.vertices.shape[0]

    @property
    def n_faces(self) -> int:
        return self.faces.shape[0]


def _normalize_to_radius(vertices: np.ndarray, radius: float) -> np.ndarray:
    norms = np.linalg.norm(vertices, axis=1)
    if np.any(norms <= np.finfo(float).tiny):
        raise ValueError("cannot normalize a mesh containing the origin")
    return radius * vertices / norms[:, None]


@dataclass(frozen=True)
class LevelSetSurface:
    """Implicit surface described by a reference mesh and level-set functions."""

    mesh: SurfaceMesh
    level_set: Callable[[np.ndarray], float]
    gradient: Callable[[np.ndarray], np.ndarray]

    @classmethod
    def sphere(
        cls,
        *,
        radius: float = 1.0,
        center: Optional[Sequence[float]] = None,
        mesh: Optional[SurfaceMesh] = None,
        mesh_refinement_level: int = 1,
    ) -> "LevelSetSurface":
        """Create a spherical level-set surface with a matching default mesh."""
        if radius <= 0.0:
            raise ValueError("radius must be positive")

        center_array = np.zeros(3, dtype=float) if center is None else np.asarray(center, dtype=float)
        if center_array.shape != (3,):
            raise ValueError("center must have shape (3,)")

        if mesh is None:
            mesh = SurfaceMesh.icosphere(
                refinement_level=mesh_refinement_level,
                radius=radius,
            )
            mesh = SurfaceMesh(mesh.vertices + center_array, mesh.faces)

        def level_set(point: np.ndarray) -> float:
            shifted = np.asarray(point, dtype=float) - center_array
            return float(np.dot(shifted, shifted) - radius**2)

        def gradient(point: np.ndarray) -> np.ndarray:
            return 2.0 * (np.asarray(point, dtype=float) - center_array)

        return cls(mesh, level_set, gradient)

    @classmethod
    def unit_sphere(
        cls,
        *,
        mesh: Optional[SurfaceMesh] = None,
        mesh_refinement_level: int = 1,
    ) -> "LevelSetSurface":
        """Create the unit sphere with a matching default icosphere mesh."""
        return cls.sphere(
            radius=1.0,
            mesh=mesh,
            mesh_refinement_level=mesh_refinement_level,
        )


@dataclass(frozen=True)
class IntegrationConfig:
    """Numerical configuration for a surface integration run."""

    interpolation_degree: int
    lp_degree: float = float("inf")
    refinement_level: int = 0
    integration_degree: int = DEFAULT_INTEGRATION_DEGREE
    quadrature_rule: str = DEFAULT_QUADRATURE_RULE

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


def _accumulate_weighted_values(
    point_values: np.ndarray,
    weights: np.ndarray,
    offsets: np.ndarray,
) -> np.ndarray:
    weighted_values = point_values * weights
    n_faces = len(offsets) - 1
    face_values = np.zeros(n_faces, dtype=float)
    if n_faces == 0:
        return face_values

    lengths = np.diff(offsets)
    if np.all(lengths > 0):
        return np.add.reduceat(weighted_values, offsets[:-1])[:n_faces]

    for face_index, (start, stop) in enumerate(zip(offsets[:-1], offsets[1:])):
        if stop > start:
            face_values[face_index] = float(np.sum(weighted_values[start:stop]))
    return face_values


@dataclass(frozen=True)
class CompiledIntegrationScheme:
    """Reusable surface quadrature for repeated integrations.

    Building the high-order geometry is usually the expensive part of a smooth
    surface integration. This object stores the quadrature points and weights
    so many scalar fields can be integrated over the same surface/configuration
    without recomputing projections and polynomial geometry.
    """

    points: np.ndarray
    weights: np.ndarray
    offsets: np.ndarray
    config: IntegrationConfig

    @property
    def n_quadrature_points(self) -> int:
        """Return the number of stored quadrature points."""
        return self.points.shape[0]

    @property
    def n_faces(self) -> int:
        """Return the number of input mesh faces represented by the scheme."""
        return len(self.offsets) - 1

    def integrate_values(self, point_values: np.ndarray) -> IntegrationResult:
        """Integrate values already evaluated at ``points``.

        Parameters
        ----------
        point_values
            One scalar value per stored quadrature point.
        """
        point_values = np.asarray(point_values, dtype=float)
        if point_values.shape != (self.n_quadrature_points,):
            raise ValueError(
                "point_values must have shape "
                f"({self.n_quadrature_points},), got {point_values.shape}"
            )

        values = _accumulate_weighted_values(point_values, self.weights, self.offsets)
        return IntegrationResult(values, self.points, self.weights, self.offsets, self.config)

    def integrate(
        self,
        integrand: Callable[[np.ndarray], float],
        *,
        vectorized: bool = False,
    ) -> IntegrationResult:
        """Integrate a scalar function using the stored quadrature.

        Set ``vectorized=True`` when ``integrand`` accepts the full point array
        with shape ``(n_points, 3)`` and returns one scalar per row.
        """
        if vectorized:
            point_values = np.asarray(integrand(self.points), dtype=float)
            point_values = np.squeeze(point_values)
        else:
            point_values = np.fromiter(
                (integrand(point) for point in self.points),
                dtype=float,
                count=self.n_quadrature_points,
            )
        return self.integrate_values(point_values)


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


@dataclass(frozen=True)
class IndicatorRefinementIteration:
    """One indicator-based reference-mesh refinement step."""

    iteration: int
    n_faces: int
    n_marked_faces: int
    max_indicator: float
    threshold: float
    marked_faces: np.ndarray
    indicator_values: Optional[np.ndarray] = None
    min_indicator: float = 0.0
    mean_indicator: float = 0.0
    median_indicator: float = 0.0
    std_indicator: float = 0.0
    marked_fraction: float = 0.0
    n_vertices: int = 0
    n_vertices_after: int = 0
    n_faces_after: int = 0
    n_new_vertices: int = 0
    n_new_faces: int = 0
    stop_reason: str = ""

    @property
    def n_unmarked_faces(self) -> int:
        """Return the number of faces below the marking threshold."""
        return self.n_faces - self.n_marked_faces

    @property
    def did_refine(self) -> bool:
        """Return whether this iteration produced a refined mesh."""
        return self.n_faces_after > self.n_faces or self.n_vertices_after > self.n_vertices


@dataclass(frozen=True)
class IndicatorRefinementResult:
    """Result returned by :func:`refine_by_indicator`."""

    final_surface: LevelSetSurface
    history: Tuple[IndicatorRefinementIteration, ...]

    @property
    def n_iterations(self) -> int:
        """Return the number of indicator refinement iterations."""
        return len(self.history)

    @property
    def n_faces(self) -> int:
        """Return the number of faces in the refined mesh."""
        return self.final_surface.mesh.n_faces

    @property
    def n_vertices(self) -> int:
        """Return the number of vertices in the refined mesh."""
        return self.final_surface.mesh.n_vertices

    @property
    def initial_n_faces(self) -> int:
        """Return the number of faces at the start of refinement."""
        if not self.history:
            return self.n_faces
        return self.history[0].n_faces

    @property
    def initial_n_vertices(self) -> int:
        """Return the number of vertices at the start of refinement."""
        if not self.history:
            return self.n_vertices
        return self.history[0].n_vertices

    @property
    def n_refined_iterations(self) -> int:
        """Return the number of iterations that actually changed the mesh."""
        return sum(1 for step in self.history if step.did_refine)

    @property
    def total_marked_faces(self) -> int:
        """Return the total number of marked faces across all iterations."""
        return sum(step.n_marked_faces for step in self.history)

    @property
    def face_growth_factor(self) -> float:
        """Return final face count divided by initial face count."""
        if self.initial_n_faces == 0:
            return 0.0
        return self.n_faces / self.initial_n_faces

    @property
    def vertex_growth_factor(self) -> float:
        """Return final vertex count divided by initial vertex count."""
        if self.initial_n_vertices == 0:
            return 0.0
        return self.n_vertices / self.initial_n_vertices

    @property
    def max_marked_fraction(self) -> float:
        """Return the largest fraction of marked faces in any iteration."""
        if not self.history:
            return 0.0
        return float(max(step.marked_fraction for step in self.history))

    @property
    def stop_reason(self) -> str:
        """Return the terminal reason reported by the last iteration."""
        if not self.history:
            return "not_started"
        return self.history[-1].stop_reason or "completed"

    @property
    def face_counts(self) -> np.ndarray:
        """Return face counts before each iteration plus the final count."""
        if not self.history:
            return np.array([self.n_faces], dtype=int)
        return np.array([step.n_faces for step in self.history] + [self.n_faces], dtype=int)

    @property
    def vertex_counts(self) -> np.ndarray:
        """Return vertex counts before each iteration plus the final count."""
        if not self.history:
            return np.array([self.n_vertices], dtype=int)
        return np.array([step.n_vertices for step in self.history] + [self.n_vertices], dtype=int)

    @property
    def max_indicators(self) -> np.ndarray:
        """Return the maximum indicator value from each iteration."""
        return np.array([step.max_indicator for step in self.history], dtype=float)

    @property
    def mean_indicators(self) -> np.ndarray:
        """Return the mean indicator value from each iteration."""
        return np.array([step.mean_indicator for step in self.history], dtype=float)

    @property
    def marked_fractions(self) -> np.ndarray:
        """Return the marked-face fraction from each iteration."""
        return np.array([step.marked_fraction for step in self.history], dtype=float)

    def summary(self) -> str:
        """Return a compact text report for the refinement run."""
        lines = [
            f"Iterations:            {self.n_iterations}",
            f"Stop reason:           {self.stop_reason}",
            f"Initial faces:         {self.initial_n_faces}",
            f"Final faces:           {self.n_faces}",
            f"Face growth factor:    {self.face_growth_factor:.3g}",
            f"Initial vertices:      {self.initial_n_vertices}",
            f"Final vertices:        {self.n_vertices}",
            f"Vertex growth factor:  {self.vertex_growth_factor:.3g}",
            f"Refined iterations:    {self.n_refined_iterations}",
            f"Total marked faces:    {self.total_marked_faces}",
            f"Max marked fraction:   {self.max_marked_fraction:.3g}",
        ]
        if self.history:
            last = self.history[-1]
            lines.extend(
                [
                    f"Last max indicator:    {last.max_indicator:.3e}",
                    f"Last mean indicator:   {last.mean_indicator:.3e}",
                    f"Last std indicator:    {last.std_indicator:.3e}",
                    f"Last threshold:        {last.threshold:.3e}",
                    f"Last marked faces:     {last.n_marked_faces}",
                    f"Last new faces:        {last.n_new_faces}",
                    f"Last new vertices:     {last.n_new_vertices}",
                ]
            )
        return "\n".join(lines)


@dataclass(frozen=True)
class PatchGeometry:
    """Geometry samples belonging to one curved surface patch."""

    index: int
    points: np.ndarray
    weights: np.ndarray
    tangent_u: np.ndarray
    tangent_v: np.ndarray
    normal: np.ndarray
    metric_tensor: np.ndarray
    area_density: np.ndarray
    second_fundamental_form: np.ndarray
    mean_curvature: np.ndarray
    gaussian_curvature: np.ndarray

    @property
    def n_points(self) -> int:
        """Return the number of sampled points on this patch."""
        return self.points.shape[0]

    @property
    def area(self) -> float:
        """Return the quadrature area of this patch."""
        return float(np.sum(self.weights))

    @property
    def center(self) -> np.ndarray:
        """Return the area-weighted center of sampled patch geometry."""
        if self.n_points == 0:
            return np.full(3, np.nan)
        area = self.area
        if area > np.finfo(float).tiny:
            return np.average(self.points, axis=0, weights=self.weights)
        return np.mean(self.points, axis=0)

    @property
    def radius(self) -> float:
        """Return the maximum sampled distance from the patch center."""
        if self.n_points == 0:
            return 0.0
        return float(np.max(np.linalg.norm(self.points - self.center, axis=1)))


@dataclass(frozen=True)
class SurfaceGeometryResult:
    """Differential geometry samples on the interpolated curved surface."""

    points: np.ndarray
    weights: np.ndarray
    offsets: np.ndarray
    tangent_u: np.ndarray
    tangent_v: np.ndarray
    normal: np.ndarray
    metric_tensor: np.ndarray
    area_density: np.ndarray
    second_fundamental_form: np.ndarray
    mean_curvature: np.ndarray
    gaussian_curvature: np.ndarray
    config: IntegrationConfig

    @property
    def n_points(self) -> int:
        """Return the number of sampled quadrature points."""
        return self.points.shape[0]

    @property
    def n_faces(self) -> int:
        """Return the number of input mesh faces represented by the offsets."""
        return len(self.offsets) - 1

    @property
    def n_patches(self) -> int:
        """Return the number of curved patches represented by the result."""
        return self.n_faces

    def patch(self, index: int) -> PatchGeometry:
        """Return geometry samples for one curved patch."""
        if index < 0 or index >= self.n_patches:
            raise IndexError("patch index out of range")

        start = self.offsets[index]
        stop = self.offsets[index + 1]
        return PatchGeometry(
            index=index,
            points=self.points[start:stop],
            weights=self.weights[start:stop],
            tangent_u=self.tangent_u[start:stop],
            tangent_v=self.tangent_v[start:stop],
            normal=self.normal[start:stop],
            metric_tensor=self.metric_tensor[start:stop],
            area_density=self.area_density[start:stop],
            second_fundamental_form=self.second_fundamental_form[start:stop],
            mean_curvature=self.mean_curvature[start:stop],
            gaussian_curvature=self.gaussian_curvature[start:stop],
        )

    @property
    def patches(self) -> Tuple[PatchGeometry, ...]:
        """Return all curved patches as patch-level geometry objects."""
        return tuple(self.patch(index) for index in range(self.n_patches))

    @property
    def patch_areas(self) -> np.ndarray:
        """Return quadrature areas for all curved patches."""
        return np.array([self.patch(index).area for index in range(self.n_patches)])

    @property
    def patch_centers(self) -> np.ndarray:
        """Return area-weighted centers for all curved patches."""
        return np.array([self.patch(index).center for index in range(self.n_patches)])

    @property
    def patch_radii(self) -> np.ndarray:
        """Return sampled patch radii around the patch centers."""
        return np.array([self.patch(index).radius for index in range(self.n_patches)])


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

    scheme = compile_integration(surface, config)
    return scheme.integrate(integrand)


def compile_integration(
    surface: LevelSetSurface,
    config: Optional[IntegrationConfig] = None,
) -> CompiledIntegrationScheme:
    """Precompute reusable quadrature for repeated integrations.

    Use this when the surface and numerical configuration are fixed but many
    scalar functions need to be integrated. The returned scheme stores the
    quadrature points, weights, offsets, and configuration.
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
        lambda _: 1.0,
        config.integration_degree,
        config.quadrature_rule,
    )
    return CompiledIntegrationScheme(points, weights, offsets, config)


def surface_geometry(
    surface: LevelSetSurface,
    config: Optional[IntegrationConfig] = None,
) -> SurfaceGeometryResult:
    """Evaluate geometry tensors and curvatures of the interpolated surface.

    The routine builds the same high-order Minterpy surface map used by
    :func:`integrate`, then differentiates the interpolant to obtain tangents,
    metric tensor, unit normal, second fundamental form, mean curvature, and
    Gaussian curvature at the quadrature points.
    """
    if config is None:
        config = IntegrationConfig(interpolation_degree=6)
    if config.interpolation_degree < 2:
        raise ValueError("interpolation_degree must be at least 2 to compute curvature")

    (
        points,
        weights,
        offsets,
        tangent_u,
        tangent_v,
        normal,
        metric_tensor,
        area_density,
        second_fundamental_form,
        mean_curvature,
        gaussian_curvature,
    ) = compute_surf_geometry(
        surface.level_set,
        surface.gradient,
        surface.mesh.vertices,
        surface.mesh.faces,
        config.interpolation_degree,
        config.lp_degree,
        config.refinement_level,
        config.integration_degree,
        config.quadrature_rule,
    )
    return SurfaceGeometryResult(
        points=points,
        weights=weights,
        offsets=offsets,
        tangent_u=tangent_u,
        tangent_v=tangent_v,
        normal=normal,
        metric_tensor=metric_tensor,
        area_density=area_density,
        second_fundamental_form=second_fundamental_form,
        mean_curvature=mean_curvature,
        gaussian_curvature=gaussian_curvature,
        config=config,
    )


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


def _face_centers(mesh: SurfaceMesh) -> np.ndarray:
    if mesh.faces.shape[1] != 3:
        raise ValueError("indicator refinement currently requires triangular faces")
    return np.mean(mesh.vertices[mesh.faces], axis=1)


def refine_by_indicator(
    surface: LevelSetSurface,
    indicator: Callable[[np.ndarray], float],
    *,
    max_iterations: int = 6,
    threshold_fraction: float = 0.25,
    use_absolute: bool = True,
) -> IndicatorRefinementResult:
    """Refine the reference mesh using an indicator evaluated at face centers.

    This reproduces the common curved-grid workflow where the linear host mesh
    is adapted first and a polynomial degree study is run afterwards on the
    adapted mesh. At each iteration, the indicator is evaluated at the affine
    center of every triangular face. Faces with indicator value larger than
    ``threshold_fraction * max_indicator`` are subdivided. Neighboring faces
    with hanging edges are split by green refinement so the adapted reference
    mesh remains conforming.
    """
    if max_iterations < 1:
        raise ValueError("max_iterations must be at least 1")
    if not (0.0 < threshold_fraction <= 1.0):
        raise ValueError("threshold_fraction must be in the interval (0, 1]")

    current_surface = surface
    history = []

    for iteration in range(max_iterations):
        centers = _face_centers(current_surface.mesh)
        indicator_values = np.asarray([indicator(center) for center in centers], dtype=float)
        if use_absolute:
            indicator_values = np.abs(indicator_values)
        if not np.all(np.isfinite(indicator_values)):
            raise ValueError("indicator returned non-finite values")

        n_faces_before = current_surface.mesh.n_faces
        n_vertices_before = current_surface.mesh.n_vertices
        if indicator_values.size:
            min_indicator = float(np.min(indicator_values))
            max_indicator = float(np.max(indicator_values))
            mean_indicator = float(np.mean(indicator_values))
            median_indicator = float(np.median(indicator_values))
            std_indicator = float(np.std(indicator_values))
        else:
            min_indicator = 0.0
            max_indicator = 0.0
            mean_indicator = 0.0
            median_indicator = 0.0
            std_indicator = 0.0
        threshold = threshold_fraction * max_indicator
        if max_indicator == 0.0:
            marked_faces = np.array([], dtype=int)
        else:
            marked_faces = np.flatnonzero(indicator_values > threshold)
        marked_fraction = marked_faces.size / n_faces_before if n_faces_before else 0.0

        if indicator_values.size == 0:
            stop_reason = "empty_mesh"
        elif max_indicator == 0.0:
            stop_reason = "zero_indicator"
        elif marked_faces.size == 0:
            stop_reason = "no_marked_faces"
        elif iteration == max_iterations - 1:
            stop_reason = "max_iterations"
        else:
            stop_reason = ""

        if stop_reason:
            history.append(
                IndicatorRefinementIteration(
                    iteration=iteration,
                    n_faces=n_faces_before,
                    n_marked_faces=int(marked_faces.size),
                    max_indicator=max_indicator,
                    threshold=float(threshold),
                    marked_faces=marked_faces,
                    indicator_values=indicator_values,
                    min_indicator=min_indicator,
                    mean_indicator=mean_indicator,
                    median_indicator=median_indicator,
                    std_indicator=std_indicator,
                    marked_fraction=float(marked_fraction),
                    n_vertices=n_vertices_before,
                    n_vertices_after=n_vertices_before,
                    n_faces_after=n_faces_before,
                    n_new_vertices=0,
                    n_new_faces=0,
                    stop_reason=stop_reason,
                )
            )
            return IndicatorRefinementResult(
                final_surface=current_surface,
                history=tuple(history),
            )

        vertices, faces = subdivide_conforming(
            current_surface.mesh.vertices,
            current_surface.mesh.faces,
            face_index=marked_faces,
        )
        current_surface = LevelSetSurface(
            SurfaceMesh(vertices, faces),
            current_surface.level_set,
            current_surface.gradient,
        )
        history.append(
            IndicatorRefinementIteration(
                iteration=iteration,
                n_faces=n_faces_before,
                n_marked_faces=int(marked_faces.size),
                max_indicator=max_indicator,
                threshold=float(threshold),
                marked_faces=marked_faces,
                indicator_values=indicator_values,
                min_indicator=min_indicator,
                mean_indicator=mean_indicator,
                median_indicator=median_indicator,
                std_indicator=std_indicator,
                marked_fraction=float(marked_fraction),
                n_vertices=n_vertices_before,
                n_vertices_after=current_surface.mesh.n_vertices,
                n_faces_after=current_surface.mesh.n_faces,
                n_new_vertices=current_surface.mesh.n_vertices - n_vertices_before,
                n_new_faces=current_surface.mesh.n_faces - n_faces_before,
            )
        )

    return IndicatorRefinementResult(
        final_surface=current_surface,
        history=tuple(history),
    )
