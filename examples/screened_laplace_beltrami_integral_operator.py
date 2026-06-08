"""Screened Laplace-Beltrami prototype using a sphere Green operator.

This example applies the exact Green operator for

    (alpha - Delta_Gamma) u = f

on the implicit unit sphere.  For the sphere, the screened Green operator is
diagonal in spherical harmonics:

    G_alpha(x, y) = sum_l (2l+1)/(4*pi*(alpha + l(l+1))) P_l(x . y).

The right-hand side is chosen as a combination of degree-one and degree-two
zonal harmonics, so the exact solution is known and the error can be measured.
This is a sphere-specific prototype, not a general surface PDE solver.
"""

from pathlib import Path

import numpy as np

from surfgeopy import (
    IntegrationConfig,
    LevelSetSurface,
    SurfaceMesh,
    surface_geometry,
)


def phi(x: np.ndarray) -> float:
    return x[0] ** 2 + x[1] ** 2 + x[2] ** 2 - 1.0


def dphi(x: np.ndarray) -> np.ndarray:
    return np.array([2.0 * x[0], 2.0 * x[1], 2.0 * x[2]])


def right_hand_side(points: np.ndarray) -> np.ndarray:
    z = points[:, 2]
    y1 = z
    y2 = 0.5 * (3.0 * z**2 - 1.0)
    return y1 + 0.25 * y2


def exact_solution(points: np.ndarray, alpha: float) -> np.ndarray:
    z = points[:, 2]
    y1 = z
    y2 = 0.5 * (3.0 * z**2 - 1.0)
    return y1 / (alpha + 2.0) + 0.25 * y2 / (alpha + 6.0)


def screened_green_kernel(
    cosine: np.ndarray,
    alpha: float,
    max_degree: int,
) -> np.ndarray:
    """Evaluate the truncated screened Green kernel on the unit sphere."""
    if alpha <= 0.0:
        raise ValueError("alpha must be positive for the screened operator")
    if max_degree < 0:
        raise ValueError("max_degree must be non-negative")

    kernel = np.zeros_like(cosine)

    p_l_minus_two = np.ones_like(cosine)
    kernel += (1.0 / (4.0 * np.pi)) / alpha * p_l_minus_two

    if max_degree == 0:
        return kernel

    p_l_minus_one = cosine
    kernel += (3.0 / (4.0 * np.pi)) / (alpha + 2.0) * p_l_minus_one

    for ell in range(2, max_degree + 1):
        p_l = (
            (2 * ell - 1) * cosine * p_l_minus_one
            - (ell - 1) * p_l_minus_two
        ) / ell
        kernel += (
            (2 * ell + 1)
            / (4.0 * np.pi)
            / (alpha + ell * (ell + 1))
            * p_l
        )
        p_l_minus_two, p_l_minus_one = p_l_minus_one, p_l

    return kernel


def apply_screened_green_operator(
    targets: np.ndarray,
    sources: np.ndarray,
    weights: np.ndarray,
    values: np.ndarray,
    alpha: float,
    max_degree: int,
    chunk_size: int = 512,
) -> np.ndarray:
    """Apply the truncated Green operator without forming a full dense matrix."""
    weighted_values = weights * values
    result = np.empty(targets.shape[0], dtype=float)

    for start in range(0, targets.shape[0], chunk_size):
        stop = min(start + chunk_size, targets.shape[0])
        cosine = np.clip(targets[start:stop] @ sources.T, -1.0, 1.0)
        kernel = screened_green_kernel(cosine, alpha, max_degree)
        result[start:stop] = kernel @ weighted_values

    return result


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    mesh = SurfaceMesh.from_mat(str(root / "meshes" / "SphereMesh_N=124_r=1.mat"))
    surface = LevelSetSurface(mesh, phi, dphi)

    config = IntegrationConfig(
        interpolation_degree=4,
        integration_degree=4,
        quadrature_rule="ModePy_VioreanuRokhlin",
    )
    geometry = surface_geometry(surface, config)
    points = geometry.points
    weights = geometry.weights

    alpha = 1.5
    rhs = right_hand_side(points)
    exact = exact_solution(points, alpha)

    print("Screened Laplace-Beltrami prototype by a Green integral operator")
    print("Equation:              (alpha - Delta_Gamma) u = f")
    print(f"alpha:                 {alpha:.3f}")
    print(f"Quadrature nodes:      {points.shape[0]}")
    print(f"Surface area check:    {np.sum(weights):.12e}")
    print("\nmax degree    relative L2 error    max error")

    for max_degree in [0, 1, 2, 3, 4, 6, 8]:
        numerical = apply_screened_green_operator(
            points,
            points,
            weights,
            rhs,
            alpha,
            max_degree,
        )
        relative_error = np.linalg.norm(numerical - exact) / np.linalg.norm(exact)
        max_error = np.max(np.abs(numerical - exact))
        print(f"{max_degree:10d}    {relative_error:17.3e}    {max_error:.3e}")

    numerical = apply_screened_green_operator(points, points, weights, rhs, alpha, 2)
    sample = np.array(
        [
            points[0],
            points[points.shape[0] // 3],
            points[2 * points.shape[0] // 3],
        ]
    )
    sample_exact = exact_solution(sample, alpha)
    sample_numerical = apply_screened_green_operator(sample, points, weights, rhs, alpha, 2)

    print("\nrepresentative surface points")
    print("point                         numerical       exact           abs error")
    for point, value, exact_value in zip(sample, sample_numerical, sample_exact):
        error = abs(value - exact_value)
        print(f"{point!s:29s} {value: .12e}  {exact_value: .12e}  {error:.3e}")

    final_error = np.linalg.norm(numerical - exact) / np.linalg.norm(exact)
    print(f"\nSelected degree-2 relative error: {final_error:.3e}")


if __name__ == "__main__":
    main()
