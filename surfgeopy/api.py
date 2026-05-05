"""High-level API for integrating functions over implicit surfaces."""

from dataclasses import dataclass
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
    "integrate",
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
