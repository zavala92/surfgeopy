"""Screened Laplace-Beltrami convergence benchmark on the unit sphere.

This benchmark is a focused prototype for the equation

    (alpha - Delta_Gamma) u = f

on the implicit unit sphere.  The manufactured solution is

    u(x, y, z) = z,

so ``Delta_Gamma u = -2u`` and ``f = (alpha + 2) u``.

The script separates three effects as far as the current prototype permits:

1. Geometry refinement: keep the Green operator and quadrature order fixed,
   vary uniform reference-mesh refinement.
2. Quadrature refinement: keep the geometry fixed, vary the quadrature order.
3. Smooth-remainder error: study the Chebyshev approximation of the smooth
   remainder in the screened Green parametrix split, independently of surface
   quadrature.

The first two studies use the exact degree-one spherical Green action.  This
avoids a dense all-pairs Green matrix: for this manufactured solution the
operator can be evaluated from the discrete first and second moments of the
quadrature rule.  The third study is a kernel-level test because the current
general-surface API does not yet expose arbitrary-surface smooth remainders.

This is a benchmark and demonstration, not a general surface PDE solver.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
from numpy.polynomial.chebyshev import chebfit, chebval
from scipy import special

from surfgeopy import (
    IntegrationConfig,
    LevelSetSurface,
    SurfaceMesh,
    surface_geometry,
)


DEFAULT_ALPHA = 0.1
DEFAULT_QUADRATURE_RULE = "ModePy_VioreanuRokhlin"


@dataclass(frozen=True)
class ErrorRow:
    """One row of a convergence study."""

    label: str
    parameter: int
    interpolation_degree: int
    refinement_level: int
    integration_degree: int
    panels: int
    points: int
    l2_error: float
    max_error: float
    rate: Optional[float] = None


@dataclass(frozen=True)
class RemainderRow:
    """One row of the smooth-remainder kernel study."""

    degree: int
    max_error: float
    max_relative_error: float


def phi(x: np.ndarray) -> float:
    return x[0] ** 2 + x[1] ** 2 + x[2] ** 2 - 1.0


def dphi(x: np.ndarray) -> np.ndarray:
    return np.array([2.0 * x[0], 2.0 * x[1], 2.0 * x[2]])


def manufactured_solution(points: np.ndarray) -> np.ndarray:
    """Return u=x_z on the sphere."""
    return np.asarray(points)[:, 2]


def manufactured_rhs(points: np.ndarray, alpha: float) -> np.ndarray:
    """Return f=(alpha+2)u for u=x_z."""
    return (alpha + 2.0) * manufactured_solution(points)


def build_surface(mesh_path: Path) -> LevelSetSurface:
    mesh = SurfaceMesh.from_mat(str(mesh_path))
    return LevelSetSurface(mesh, phi, dphi)


def geometry_for_config(
    surface: LevelSetSurface,
    interpolation_degree: int,
    refinement_level: int,
    integration_degree: int,
    quadrature_rule: str = DEFAULT_QUADRATURE_RULE,
):
    config = IntegrationConfig(
        interpolation_degree=interpolation_degree,
        refinement_level=refinement_level,
        integration_degree=integration_degree,
        quadrature_rule=quadrature_rule,
    )
    return surface_geometry(surface, config)


def apply_degree_one_screened_green(
    points: np.ndarray,
    weights: np.ndarray,
    alpha: float,
) -> np.ndarray:
    """Apply the exact degree-one screened sphere Green action to f=(alpha+2)z."""
    rhs_z_moment = float(np.dot(weights, points[:, 2]))
    yz_moment = points.T @ (weights * points[:, 2])

    constant_part = ((alpha + 2.0) / (4.0 * np.pi * alpha)) * rhs_z_moment
    degree_one_part = (3.0 / (4.0 * np.pi)) * (points @ yz_moment)
    return constant_part + degree_one_part


def weighted_errors(points: np.ndarray, weights: np.ndarray, numerical: np.ndarray) -> tuple[float, float]:
    exact = manufactured_solution(points)
    error = numerical - exact
    denominator = np.sqrt(np.dot(weights, exact * exact))
    l2_error = np.sqrt(np.dot(weights, error * error)) / max(denominator, np.finfo(float).tiny)
    max_error = float(np.max(np.abs(error)))
    return float(l2_error), max_error


def effective_panels(surface: LevelSetSurface, refinement_level: int) -> int:
    return int(surface.mesh.n_faces * (4 ** refinement_level))


def with_rates(rows: list[ErrorRow]) -> list[ErrorRow]:
    rated = []
    previous = None
    for row in rows:
        rate = None
        if previous is not None and row.l2_error > 0.0 and previous.l2_error > 0.0:
            h_previous = 1.0 / np.sqrt(previous.panels)
            h_current = 1.0 / np.sqrt(row.panels)
            if h_current != h_previous:
                rate = np.log(previous.l2_error / row.l2_error) / np.log(
                    h_previous / h_current
                )
        rated.append(ErrorRow(**{**row.__dict__, "rate": None if rate is None else float(rate)}))
        previous = row
    return rated


def run_geometry_study(
    mesh_path: Path,
    alpha: float = DEFAULT_ALPHA,
    interpolation_degree: int = 4,
    integration_degree: int = 6,
    refinement_levels: Iterable[int] = (0, 1, 2),
) -> list[ErrorRow]:
    surface = build_surface(mesh_path)
    rows = []
    for refinement_level in refinement_levels:
        geometry = geometry_for_config(
            surface,
            interpolation_degree,
            refinement_level,
            integration_degree,
        )
        numerical = apply_degree_one_screened_green(geometry.points, geometry.weights, alpha)
        l2_error, max_error = weighted_errors(geometry.points, geometry.weights, numerical)
        rows.append(
            ErrorRow(
                label="geometry",
                parameter=int(refinement_level),
                interpolation_degree=int(interpolation_degree),
                refinement_level=int(refinement_level),
                integration_degree=int(integration_degree),
                panels=effective_panels(surface, refinement_level),
                points=geometry.n_points,
                l2_error=l2_error,
                max_error=max_error,
            )
        )
    return with_rates(rows)


def run_quadrature_study(
    mesh_path: Path,
    alpha: float = DEFAULT_ALPHA,
    interpolation_degree: int = 4,
    refinement_level: int = 1,
    integration_degrees: Iterable[int] = (1, 2, 3, 4, 6),
) -> list[ErrorRow]:
    surface = build_surface(mesh_path)
    rows = []
    for integration_degree in integration_degrees:
        geometry = geometry_for_config(
            surface,
            interpolation_degree,
            refinement_level,
            integration_degree,
        )
        numerical = apply_degree_one_screened_green(geometry.points, geometry.weights, alpha)
        l2_error, max_error = weighted_errors(geometry.points, geometry.weights, numerical)
        rows.append(
            ErrorRow(
                label="quadrature",
                parameter=int(integration_degree),
                interpolation_degree=int(interpolation_degree),
                refinement_level=int(refinement_level),
                integration_degree=int(integration_degree),
                panels=effective_panels(surface, refinement_level),
                points=geometry.n_points,
                l2_error=l2_error,
                max_error=max_error,
            )
        )
    return rows


def exact_sphere_green(cosine: np.ndarray, alpha: float) -> np.ndarray:
    """Exact screened Green kernel on the unit sphere for 0<alpha<1/4."""
    if alpha <= 0.0 or alpha >= 0.25:
        raise ValueError("the closed-form kernel study requires 0 < alpha < 1/4")
    degree = (-1.0 + np.sqrt(1.0 - 4.0 * alpha)) / 2.0
    return -special.lpmv(0, degree, -cosine) / (4.0 * np.sin(np.pi * degree))


def curvature_corrected_parametrix(cosine: np.ndarray, alpha: float) -> np.ndarray:
    chord = np.sqrt(np.maximum(2.0 * (1.0 - cosine), np.finfo(float).tiny))
    distance = chord * (1.0 + chord**2 / 24.0 + 3.0 * chord**4 / 640.0)
    return special.k0(np.sqrt(alpha) * distance) / (2.0 * np.pi)


def fit_smooth_remainder(alpha: float, degree: int, n_samples: int = 1600) -> np.ndarray:
    theta = np.linspace(np.pi, 1.0e-4, n_samples)
    cosine = np.cos(theta)
    exact = exact_sphere_green(cosine, alpha)
    singular = curvature_corrected_parametrix(cosine, alpha)
    return chebfit(cosine, exact - singular, degree)


def run_remainder_study(
    alpha: float = DEFAULT_ALPHA,
    degrees: Iterable[int] = (2, 4, 8, 12, 20, 30),
) -> list[RemainderRow]:
    theta = np.geomspace(1.0e-3, np.pi, 1000)
    cosine = np.cos(theta)
    exact = exact_sphere_green(cosine, alpha)
    singular = curvature_corrected_parametrix(cosine, alpha)
    rows = []
    for degree in degrees:
        coefficients = fit_smooth_remainder(alpha, degree)
        reconstructed = singular + chebval(cosine, coefficients)
        error = np.abs(reconstructed - exact)
        rows.append(
            RemainderRow(
                degree=int(degree),
                max_error=float(np.max(error)),
                max_relative_error=float(np.max(error / np.maximum(np.abs(exact), 1.0e-30))),
            )
        )
    return rows


def print_error_table(title: str, rows: list[ErrorRow], parameter_name: str) -> None:
    print(f"\n{title}")
    print(
        f"{parameter_name:>8s}  {'k':>3s}  {'q':>3s}  {'panels':>8s}  "
        f"{'points':>8s}  {'rel L2':>12s}  {'max':>12s}  {'rate':>8s}"
    )
    for row in rows:
        rate = "" if row.rate is None else f"{row.rate:8.2f}"
        print(
            f"{row.parameter:8d}  "
            f"{row.interpolation_degree:3d}  "
            f"{row.integration_degree:3d}  "
            f"{row.panels:8d}  "
            f"{row.points:8d}  "
            f"{row.l2_error:12.3e}  "
            f"{row.max_error:12.3e}  "
            f"{rate:>8s}"
        )


def print_remainder_table(rows: list[RemainderRow]) -> None:
    print("\nSmooth-remainder / screened-parametrix kernel study")
    print(f"{'degree':>8s}  {'max kernel err':>16s}  {'max relative err':>16s}")
    for row in rows:
        print(
            f"{row.degree:8d}  "
            f"{row.max_error:16.3e}  "
            f"{row.max_relative_error:16.3e}"
        )


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    mesh_path = root / "meshes" / "SphereMesh_N=124_r=1.mat"
    alpha = DEFAULT_ALPHA

    print("Screened Laplace-Beltrami sphere convergence benchmark")
    print("Manufactured solution: u(x,y,z)=z")
    print("Equation:              (alpha - Delta_Gamma) u = (alpha + 2) u")
    print(f"alpha:                 {alpha:.3f}")
    print(f"mesh:                  {mesh_path}")

    geometry_rows = run_geometry_study(mesh_path, alpha)
    print_error_table(
        "A. Geometry refinement study",
        geometry_rows,
        "ref",
    )

    quadrature_rows = run_quadrature_study(mesh_path, alpha)
    print_error_table(
        "B. Quadrature refinement study",
        quadrature_rows,
        "qdeg",
    )

    remainder_rows = run_remainder_study(alpha)
    print_remainder_table(remainder_rows)

    print(
        "\nInterpretation: the first two tables measure the exact degree-one "
        "sphere Green action with surfgeopy geometry and quadrature. The third "
        "table isolates the kernel-level smooth-remainder approximation used by "
        "the screened parametrix prototype."
    )


if __name__ == "__main__":
    main()
