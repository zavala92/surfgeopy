"""Nonconstant-density error check for the prototype SSPI quadrature.

The density sigma(y) = y_z is a degree-1 spherical harmonic on the unit sphere.
For the Laplace single-layer potential with kernel 1/(4*pi*|x-y|), the exact
solution is

    u(x) = x_z / 3                 for |x| = 1,
    u(x) = direction_z / (3*r^2)   for x = r*direction, r > 1.

This tests whether the prototype singular quadrature works beyond constant
densities.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from surfgeopy import (
    IntegrationConfig,
    LaplaceSingleLayerOperator,
    LevelSetSurface,
    SingularIntegrationConfig,
    SurfaceMesh,
    integrate,
)


def phi(x: np.ndarray) -> float:
    return x[0] ** 2 + x[1] ** 2 + x[2] ** 2 - 1.0


def dphi(x: np.ndarray) -> np.ndarray:
    return np.array([2.0 * x[0], 2.0 * x[1], 2.0 * x[2]])


def density(point: np.ndarray) -> float:
    return float(point[2])


def exact_degree_one_single_layer(targets: np.ndarray) -> np.ndarray:
    radii = np.linalg.norm(targets, axis=1)
    directions = targets / radii[:, None]
    exact = directions[:, 2] / (3.0 * radii**2)
    on_or_inside = radii <= 1.0
    exact[on_or_inside] = targets[on_or_inside, 2] / 3.0
    return exact


def standard_single_layer(surface: LevelSetSurface, targets: np.ndarray) -> np.ndarray:
    config = IntegrationConfig(
        interpolation_degree=6,
        integration_degree=16,
        refinement_level=0,
        quadrature_rule="ModePy_XiaoGimbutas",
    )
    values = []
    for target in targets:
        def kernel(point: np.ndarray, target=target) -> float:
            return density(point) / (4.0 * np.pi * np.linalg.norm(point - target))

        values.append(integrate(surface, kernel, config).total)

    return np.array(values)


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    mesh = SurfaceMesh.from_mat(str(root / "meshes" / "SphereMesh_N=124_r=1.mat"))
    surface = LevelSetSurface(mesh, phi, dphi)

    direction = mesh.vertices[np.argmax(np.abs(mesh.vertices[:, 2]))]
    direction = direction / np.linalg.norm(direction)
    radii = np.array([1.0, 1.001, 1.005, 1.01, 1.05, 1.10, 1.50, 2.00])
    targets = radii[:, None] * direction[None, :]
    exact = exact_degree_one_single_layer(targets)

    sspi_config = SingularIntegrationConfig(
        interpolation_degree=6,
        regular_order=12,
        smooth_degree=8,
        moment_order=44,
        correction_order=14,
        near_threshold=1.5,
        singular_model="curvature",
    )

    standard_values = standard_single_layer(surface, targets)
    sspi_operator = LaplaceSingleLayerOperator(surface, sspi_config)
    sspi_result = sspi_operator.evaluate(targets, density=density)

    standard_error = np.abs(standard_values - exact)
    sspi_error = np.abs(sspi_result.values - exact)
    distance_to_surface = np.maximum(radii - 1.0, 1.0e-16)

    print(
        "radius   distance    standard error    SSPI error       "
        "near panels   singular panels"
    )
    for index, radius in enumerate(radii):
        print(
            f"{radius:6.3f}  "
            f"{distance_to_surface[index]:9.1e}  "
            f"{standard_error[index]:14.6e}  "
            f"{sspi_error[index]:14.6e}  "
            f"{sspi_result.near_panel_counts[index]:11d}  "
            f"{sspi_result.singular_panel_counts[index]:15d}"
        )

    print(f"\nmax standard error: {np.max(standard_error):.6e}")
    print(f"max SSPI error:     {np.max(sspi_error):.6e}")

    figure_path = root / "images" / "laplace_single_layer_spherical_harmonic_error.png"
    figure_path.parent.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(7.0, 4.8))
    plt.loglog(distance_to_surface, standard_error, "-o", label="standard quadrature")
    plt.loglog(distance_to_surface, sspi_error, "-s", label="SSPI prototype")
    plt.xlabel("Distance to unit sphere")
    plt.ylabel("Absolute error")
    plt.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(figure_path, dpi=200)
    print(f"\nSaved plot to {figure_path}")


if __name__ == "__main__":
    main()
