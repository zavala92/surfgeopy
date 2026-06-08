"""Prototype harmonic extension from an implicit surface.

This example demonstrates a small dense boundary-integral calculation for the
interior/exterior Dirichlet problem on the unit sphere with boundary data
g(x) = x_z.  The single-layer ansatz is

    u(x) = integral_Gamma sigma(y) / (4*pi*|x-y|) dS_y.

On the unit sphere the exact density is sigma(y) = 3*y_z, and the exact
harmonic extension is u(x)=x_z inside the sphere and u(x)=x_z/|x|^3 outside.
This is a validation prototype, not a reusable integral-equation solver API.
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


def boundary_data(points: np.ndarray) -> np.ndarray:
    return points[:, 2]


def exact_density(points: np.ndarray) -> np.ndarray:
    return 3.0 * points[:, 2]


def exact_harmonic_extension(targets: np.ndarray) -> np.ndarray:
    radii = np.linalg.norm(targets, axis=1)
    values = targets[:, 2] / np.maximum(radii, np.finfo(float).tiny) ** 3
    values[radii <= 1.0] = targets[radii <= 1.0, 2]
    return values


def assemble_single_layer_matrix(points: np.ndarray, weights: np.ndarray) -> np.ndarray:
    distances = np.linalg.norm(points[:, None, :] - points[None, :, :], axis=2)
    matrix = weights[None, :] / (4.0 * np.pi * np.maximum(distances, np.finfo(float).tiny))

    # Disk self-panel model: integral_{disk(a)} 1/(4*pi*r) dA = a/2.
    self_panel = 0.5 * np.sqrt(weights / np.pi)
    np.fill_diagonal(matrix, self_panel)
    return matrix


def evaluate_single_layer(
    targets: np.ndarray,
    points: np.ndarray,
    weights: np.ndarray,
    density: np.ndarray,
) -> np.ndarray:
    distances = np.linalg.norm(targets[:, None, :] - points[None, :, :], axis=2)
    kernel = 1.0 / (4.0 * np.pi * np.maximum(distances, np.finfo(float).tiny))
    return kernel @ (weights * density)


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    mesh = SurfaceMesh.from_mat(str(root / "meshes" / "SphereMesh_N=124_r=1.mat"))
    surface = LevelSetSurface(mesh, phi, dphi)

    config = IntegrationConfig(
        interpolation_degree=4,
        integration_degree=2,
        quadrature_rule="ModePy_VioreanuRokhlin",
    )
    geometry = surface_geometry(surface, config)
    points = geometry.points
    weights = geometry.weights

    matrix = assemble_single_layer_matrix(points, weights)
    rhs = boundary_data(points)
    density = np.linalg.solve(matrix, rhs)

    boundary_residual = np.linalg.norm(matrix @ density - rhs) / np.linalg.norm(rhs)
    density_error = np.linalg.norm(density - exact_density(points)) / np.linalg.norm(
        exact_density(points)
    )

    targets = np.array(
        [
            [0.2, 0.1, 0.3],
            [0.0, 0.0, 0.5],
            [0.3, -0.2, 1.4],
            [0.0, 0.0, 2.0],
        ],
        dtype=float,
    )
    numerical = evaluate_single_layer(targets, points, weights, density)
    exact = exact_harmonic_extension(targets)
    errors = np.abs(numerical - exact)

    print("Prototype harmonic extension by a single-layer boundary integral equation")
    print(f"Quadrature nodes:       {points.shape[0]}")
    print(f"Matrix condition no.:   {np.linalg.cond(matrix):.3e}")
    print(f"Boundary residual:      {boundary_residual:.3e}")
    print(f"Density relative error: {density_error:.3e}")
    print("\ntarget                  numerical       exact           abs error")
    for target, value, exact_value, error in zip(targets, numerical, exact, errors):
        print(
            f"{target!s:24s} "
            f"{value: .12e}  "
            f"{exact_value: .12e}  "
            f"{error:.3e}"
        )

    print(f"\nMax field error:        {np.max(errors):.3e}")


if __name__ == "__main__":
    main()
