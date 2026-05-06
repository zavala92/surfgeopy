"""Reference-domain quadrature helpers."""

from dataclasses import dataclass
from math import ceil
from typing import Tuple

import numpy as np

from .quadrature_points import quadrule_on_simplex
from .quadrature_points_gl import gauss_legendre_square
from .utils import pullback

__all__ = [
    "PULL_BACK_GAUSS",
    "GAUSS_LEGENDRE",
    "MODEPY_XIAO_GIMBUTAS",
    "MODEPY_GRUNDMANN_MOELLER",
    "MODEPY_VIOREANU_ROKHLIN",
    "MODEPY_SIMPLEX_RULES",
    "DEFAULT_QUADRATURE_RULE",
    "ReferenceQuadrature",
    "make_reference_quadrature",
]

PULL_BACK_GAUSS = "Pull_back_Gauss"
GAUSS_LEGENDRE = "Gauss_Legendre"
MODEPY_XIAO_GIMBUTAS = "ModePy_XiaoGimbutas"
MODEPY_GRUNDMANN_MOELLER = "ModePy_GrundmannMoeller"
MODEPY_VIOREANU_ROKHLIN = "ModePy_VioreanuRokhlin"
MODEPY_SIMPLEX_RULES = (
    MODEPY_XIAO_GIMBUTAS,
    MODEPY_GRUNDMANN_MOELLER,
    MODEPY_VIOREANU_ROKHLIN,
)
DEFAULT_QUADRATURE_RULE = MODEPY_VIOREANU_ROKHLIN


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
        if self.rule != PULL_BACK_GAUSS and self.rule not in MODEPY_SIMPLEX_RULES:
            return 1.0

        point = self.reference_points[index]
        return 8.0 / np.sqrt(
            (point[0] - point[1]) ** 2 + 4.0 * (1.0 - point[0] - point[1])
        )


def make_reference_quadrature(
    degree: int,
    rule: str = DEFAULT_QUADRATURE_RULE,
) -> ReferenceQuadrature:
    """Build quadrature data for the requested reference rule."""
    if rule == PULL_BACK_GAUSS:
        weights, reference_points = quadrule_on_simplex(degree)
        evaluation_points = pullback(reference_points, duffy_transform=False)
        return ReferenceQuadrature(weights, reference_points, evaluation_points, rule)

    if rule == GAUSS_LEGENDRE:
        weights, reference_points = gauss_legendre_square(degree)
        return ReferenceQuadrature(weights, reference_points, reference_points, rule)

    if rule in MODEPY_SIMPLEX_RULES:
        weights, reference_points = _modepy_simplex_quadrature(degree, rule)
        evaluation_points = pullback(reference_points, duffy_transform=False)
        return ReferenceQuadrature(weights, reference_points, evaluation_points, rule)

    raise ValueError(f"Unknown quadrature rule: {rule}")


def _modepy_simplex_quadrature(degree: int, rule: str) -> Tuple[np.ndarray, np.ndarray]:
    """Return a ModePy triangle rule mapped to surfgeopy's unit simplex."""
    try:
        import modepy as mp
    except ImportError as exc:
        raise ImportError(
            f"{rule} requires the modepy package. Install surfgeopy with modepy "
            "available in the active environment."
        ) from exc

    try:
        if rule == MODEPY_XIAO_GIMBUTAS:
            quadrature = mp.XiaoGimbutasSimplexQuadrature(degree, 2)
        elif rule == MODEPY_GRUNDMANN_MOELLER:
            order = max(0, ceil((degree - 1) / 2))
            quadrature = mp.GrundmannMoellerSimplexQuadrature(order, 2)
        elif rule == MODEPY_VIOREANU_ROKHLIN:
            quadrature = mp.VioreanuRokhlinSimplexQuadrature(degree, 2)
        else:
            raise ValueError(f"Unknown ModePy simplex rule: {rule}")
    except mp.QuadratureRuleUnavailable as exc:
        raise ValueError(
            f"{rule} is unavailable in ModePy for degree {degree}. "
            f"Try {MODEPY_XIAO_GIMBUTAS} for higher-degree simplex quadrature."
        ) from exc

    nodes = np.asarray(quadrature.nodes, dtype=float)
    if nodes.ndim != 2 or nodes.shape[0] != 2:
        raise ValueError(f"{rule} returned nodes with unexpected shape {nodes.shape}")

    # ModePy triangle rules use the biunit simplex. Map to surfgeopy's
    # reference triangle with vertices (0, 0), (1, 0), (0, 1).
    reference_points = (nodes.T + 1.0) / 2.0
    weights = np.asarray(quadrature.weights, dtype=float) / 4.0
    return weights, reference_points
