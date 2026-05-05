"""Reference-domain quadrature helpers."""

from dataclasses import dataclass
from typing import Tuple

import numpy as np

from .quadrature_points import quadrule_on_simplex
from .quadrature_points_gl import gauss_legendre_square
from .utils import pullback

__all__ = [
    "PULL_BACK_GAUSS",
    "GAUSS_LEGENDRE",
    "RECURSIVE_NODES_GAUSS_LEGENDRE",
    "ReferenceQuadrature",
    "make_reference_quadrature",
]

PULL_BACK_GAUSS = "Pull_back_Gauss"
GAUSS_LEGENDRE = "Gauss_Legendre"
RECURSIVE_NODES_GAUSS_LEGENDRE = "RecursiveNodes_GaussLegendre"


@dataclass(frozen=True)
class ReferenceQuadrature:
    """Quadrature data on the polynomial evaluation domain."""

    weights: np.ndarray
    reference_points: np.ndarray
    evaluation_points: np.ndarray
    rule: str

    @property
    def size(self) -> int:
        return self.weights.shape[0]

    def weight_scale(self, index: int) -> float:
        """Return the extra square-squeezing weight scale for one point."""
        if self.rule not in (PULL_BACK_GAUSS, RECURSIVE_NODES_GAUSS_LEGENDRE):
            return 1.0

        point = self.reference_points[index]
        return 8.0 / np.sqrt(
            (point[0] - point[1]) ** 2 + 4.0 * (1.0 - point[0] - point[1])
        )


def make_reference_quadrature(degree: int, rule: str = PULL_BACK_GAUSS) -> ReferenceQuadrature:
    """Build quadrature data for the requested reference rule."""
    if rule == PULL_BACK_GAUSS:
        weights, reference_points = quadrule_on_simplex(degree)
        evaluation_points = pullback(reference_points, duffy_transform=False)
        return ReferenceQuadrature(weights, reference_points, evaluation_points, rule)

    if rule == GAUSS_LEGENDRE:
        weights, reference_points = gauss_legendre_square(degree)
        return ReferenceQuadrature(weights, reference_points, reference_points, rule)

    if rule == RECURSIVE_NODES_GAUSS_LEGENDRE:
        weights, reference_points = _recursivenodes_simplex_gauss_legendre(degree)
        evaluation_points = pullback(reference_points, duffy_transform=False)
        return ReferenceQuadrature(weights, reference_points, evaluation_points, rule)

    raise ValueError(f"Unknown quadrature rule: {rule}")


def _recursivenodes_simplex_gauss_legendre(degree: int) -> Tuple[np.ndarray, np.ndarray]:
    """Return recursivenodes Gaussian quadrature on the unit triangle."""
    try:
        from recursivenodes.quadrature import simplexgausslegendre
    except ImportError as exc:
        raise ImportError(
            "RecursiveNodes_GaussLegendre requires the recursivenodes package."
        ) from exc

    biunit_points, biunit_weights = simplexgausslegendre(2, degree)
    reference_points = (biunit_points.reshape(-1, 2) + 1.0) / 2.0
    weights = biunit_weights.reshape(-1) / 4.0
    return weights, reference_points
