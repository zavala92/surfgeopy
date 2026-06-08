"""Surface geometry helpers for projection-based quadrature."""

from dataclasses import dataclass
from typing import Callable

import numpy as np

__all__ = [
    "ImplicitSurface",
    "ProjectionResult",
    "simplex_barycentric_coordinates",
    "affine_triangle_points",
    "project_triangle_nodes",
]


@dataclass(frozen=True)
class ProjectionResult:
    """Diagnostic information for one closest-point projection."""

    point: np.ndarray
    converged: bool
    iterations: int
    residual: float


@dataclass(frozen=True)
class ImplicitSurface:
    """Implicit surface represented by a level-set function and its gradient."""

    level_set: Callable[[np.ndarray], float]
    gradient: Callable[[np.ndarray], np.ndarray]
    projection_max_iter: int = 10
    projection_tolerance: float = 10 * np.finfo(np.float64).eps

    def project(self, point: np.ndarray) -> np.ndarray:
        """Project a point to the zero level set without mutating the input."""
        return self.project_with_info(point).point

    def project_with_info(self, point: np.ndarray) -> ProjectionResult:
        """Project a point and return convergence diagnostics."""
        projected = np.array(point, dtype=float, copy=True)
        residual = abs(float(self.level_set(projected)))

        for iteration in range(self.projection_max_iter + 1):
            if residual <= self.projection_tolerance:
                return ProjectionResult(projected, True, iteration, residual)
            if iteration == self.projection_max_iter:
                break

            grad_phi = np.asarray(self.gradient(projected), dtype=float)
            grad_phi_norm = float(np.dot(grad_phi, grad_phi))
            if grad_phi_norm <= np.finfo(np.float64).tiny:
                break

            projected -= (float(self.level_set(projected)) / grad_phi_norm) * grad_phi
            residual = abs(float(self.level_set(projected)))

        return ProjectionResult(projected, False, self.projection_max_iter, residual)


def simplex_barycentric_coordinates(points: np.ndarray) -> np.ndarray:
    """Convert reference-simplex coordinates `(x, y)` to barycentric weights."""
    points = np.asarray(points)
    return np.column_stack((1.0 - points[:, 0] - points[:, 1], points[:, 0], points[:, 1]))


def affine_triangle_points(vertices: np.ndarray, barycentric_points: np.ndarray) -> np.ndarray:
    """Evaluate the affine reference-triangle map at barycentric points."""
    return np.asarray(barycentric_points) @ np.asarray(vertices)


def project_triangle_nodes(
    surface: ImplicitSurface,
    vertices: np.ndarray,
    barycentric_points: np.ndarray,
) -> np.ndarray:
    """Project affine triangle nodes onto an implicit surface."""
    affine_points = affine_triangle_points(vertices, barycentric_points)
    projected_points = []
    for index, point in enumerate(affine_points):
        projection = surface.project_with_info(point)
        if not projection.converged:
            raise RuntimeError(
                "Projection failed to converge for triangle node "
                f"{index}: residual={projection.residual:.3e}, "
                f"iterations={projection.iterations}"
            )
        projected_points.append(projection.point)
    return np.array(projected_points)
