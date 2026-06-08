"""Screened Laplace-Beltrami prototype with SSPI parametrix quadrature.

This example uses the package-level screened parametrix operator.  The kernel
is split into a local singular part handled by SSPI and a smooth Chebyshev
remainder:

    G_alpha(x, y) = K0(sqrt(alpha) d(x, y))/(2*pi) + R_alpha(x, y).

The validation surface is the implicit unit sphere, where the exact screened
Green kernel and the exact solution are known.  The same machinery is an early
prototype path for general implicit manifolds once a local/global remainder
model is supplied.
"""

from pathlib import Path

import numpy as np
from numpy.polynomial.chebyshev import chebfit, chebval
from scipy import special

from surfgeopy import (
    LevelSetSurface,
    ScreenedLaplaceBeltramiParametrixOperator,
    ScreenedParametrixConfig,
    SurfaceMesh,
)


def phi(x: np.ndarray) -> float:
    return x[0] ** 2 + x[1] ** 2 + x[2] ** 2 - 1.0


def dphi(x: np.ndarray) -> np.ndarray:
    return np.array([2.0 * x[0], 2.0 * x[1], 2.0 * x[2]])


def right_hand_side(point: np.ndarray) -> float:
    z = point[2]
    y1 = z
    y2 = 0.5 * (3.0 * z**2 - 1.0)
    return y1 + 0.25 * y2


def exact_solution(points: np.ndarray, alpha: float) -> np.ndarray:
    z = points[:, 2]
    y1 = z
    y2 = 0.5 * (3.0 * z**2 - 1.0)
    return y1 / (alpha + 2.0) + 0.25 * y2 / (alpha + 6.0)


def exact_sphere_green(cosine: np.ndarray, alpha: float) -> np.ndarray:
    if alpha <= 0.0 or alpha >= 0.25:
        raise ValueError("this closed form requires 0 < alpha < 1/4")

    degree = (-1.0 + np.sqrt(1.0 - 4.0 * alpha)) / 2.0
    return -special.lpmv(0, degree, -cosine) / (4.0 * np.sin(np.pi * degree))


def singular_parametrix(cosine: np.ndarray, alpha: float) -> np.ndarray:
    chord = np.sqrt(np.maximum(2.0 * (1.0 - cosine), np.finfo(float).tiny))
    distance = chord * (1.0 + chord**2 / 24.0 + 3.0 * chord**4 / 640.0)
    return special.k0(np.sqrt(alpha) * distance) / (2.0 * np.pi)


def build_smooth_remainder(alpha: float, degree: int = 20):
    samples = np.cos(np.linspace(np.pi, 1.0e-4, 2000))
    coefficients = chebfit(
        samples,
        exact_sphere_green(samples, alpha) - singular_parametrix(samples, alpha),
        degree,
    )

    def smooth_remainder(source: np.ndarray, target: np.ndarray) -> float:
        cosine = np.clip(np.dot(source, target), -1.0, 1.0)
        return float(chebval(cosine, coefficients))

    return smooth_remainder


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    mesh = SurfaceMesh.from_mat(str(root / "meshes" / "SphereMesh_N=124_r=1.mat"))
    surface = LevelSetSurface(mesh, phi, dphi)

    alpha = 0.1
    smooth_remainder = build_smooth_remainder(alpha)
    config = ScreenedParametrixConfig(
        interpolation_degree=4,
        regular_order=10,
        smooth_degree=5,
        moment_order=18,
        near_threshold=1.5,
    )
    operator = ScreenedLaplaceBeltramiParametrixOperator(
        surface,
        alpha,
        config,
        smooth_remainder,
    )

    target_ids = np.array([0, mesh.n_vertices // 5, 2 * mesh.n_vertices // 5])
    targets = mesh.vertices[target_ids]
    targets = targets / np.linalg.norm(targets, axis=1)[:, None]

    result = operator.evaluate(targets, right_hand_side)
    exact = exact_solution(targets, alpha)
    errors = np.abs(result.values - exact)
    relative_error = np.linalg.norm(result.values - exact) / np.linalg.norm(exact)

    print("Screened Laplace-Beltrami prototype by SSPI parametrix quadrature")
    print("Equation:              (alpha - Delta_Gamma) u = f")
    print(f"alpha:                 {alpha:.3f}")
    print(f"Curved patches:        {operator.n_patches}")
    print(f"Targets:               {targets.shape[0]}")
    print(f"Relative error:        {relative_error:.3e}")
    print(f"Max abs error:         {np.max(errors):.3e}")
    print("\ntarget                  numerical       exact           abs error")

    for target, numerical, exact_value, error, near_count, tail in zip(
        targets,
        result.values,
        exact,
        errors,
        result.near_panel_counts,
        result.tail_indicators,
    ):
        print(
            f"{target!s:24s} "
            f"{numerical: .12e}  "
            f"{exact_value: .12e}  "
            f"{error:.3e}  "
            f"near={near_count:3d}  "
            f"tail={tail:.3e}"
        )


if __name__ == "__main__":
    main()
