"""Near-singular convergence check for the prototype SSPI quadrature.

The application is the Laplace single-layer potential on the unit sphere with
constant density,

    u(x) = integral_{S^2} 1 / (4*pi*|x-y|) dS_y.

The exact answer is

    u(x) = 1       for |x| <= 1,
    u(x) = 1/|x|   for |x| > 1.

This script compares ordinary high-order surface quadrature with the prototype
square-squeezed product integration (SSPI) as the target approaches the
surface.
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


def exact_single_layer_unit_sphere(targets: np.ndarray) -> np.ndarray:
    radii = np.linalg.norm(targets, axis=1)
    return np.where(radii <= 1.0, 1.0, 1.0 / radii)


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
            distance = np.linalg.norm(point - target)
            return 1.0 / (4.0 * np.pi * distance)

        values.append(integrate(surface, kernel, config).total)

    return np.array(values)


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    mesh = SurfaceMesh.from_mat(str(root / "meshes" / "SphereMesh_N=124_r=1.mat"))
    surface = LevelSetSurface(mesh, phi, dphi)

    direction = mesh.vertices[0] / np.linalg.norm(mesh.vertices[0])
    radii = np.array([1.0, 1.001, 1.005, 1.01, 1.05, 1.10, 1.50, 2.00])
    targets = radii[:, None] * direction[None, :]
    exact = exact_single_layer_unit_sphere(targets)

    sspi_config = SingularIntegrationConfig(
        interpolation_degree=6,
        regular_order=12,
        smooth_degree=6,
        moment_order=28,
        correction_order=14,
        near_threshold=1.5,
        singular_model="curvature",
    )
    sspi_accuracy_config = SingularIntegrationConfig(
        interpolation_degree=10,
        regular_order=14,
        smooth_degree=10,
        moment_order=60,
        correction_order=22,
        near_threshold=1.5,
        singular_model="curvature",
    )

    standard_values = standard_single_layer(surface, targets)
    sspi_operator = LaplaceSingleLayerOperator(surface, sspi_config)
    sspi_result = sspi_operator.evaluate(targets)
    sspi_accuracy_operator = LaplaceSingleLayerOperator(surface, sspi_accuracy_config)
    sspi_accuracy_result = sspi_accuracy_operator.evaluate(targets)

    standard_error = np.abs(standard_values - exact)
    sspi_error = np.abs(sspi_result.values - exact)
    sspi_accuracy_error = np.abs(sspi_accuracy_result.values - exact)
    distance_to_surface = np.maximum(radii - 1.0, 1.0e-16)

    print(
        "radius   distance    standard error    SSPI error       SSPI accuracy    "
        "near panels   singular panels"
    )
    for index, radius in enumerate(radii):
        print(
            f"{radius:6.3f}  "
            f"{distance_to_surface[index]:9.1e}  "
            f"{standard_error[index]:14.6e}  "
            f"{sspi_error[index]:14.6e}  "
            f"{sspi_accuracy_error[index]:14.6e}  "
            f"{sspi_result.near_panel_counts[index]:11d}  "
            f"{sspi_result.singular_panel_counts[index]:15d}"
        )

    print(f"\nmax standard error: {np.max(standard_error):.6e}")
    print(f"max SSPI error:     {np.max(sspi_error):.6e}")
    print(f"max SSPI accuracy:  {np.max(sspi_accuracy_error):.6e}")

    figure_path = root / "images" / "laplace_single_layer_near_singular_convergence.png"
    figure_path.parent.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(7.0, 4.8))
    plt.loglog(distance_to_surface, standard_error, "-o", label="standard quadrature")
    plt.loglog(distance_to_surface, sspi_error, "-s", label="SSPI prototype")
    plt.loglog(distance_to_surface, sspi_accuracy_error, "-^", label="SSPI accuracy")
    plt.xlabel("Distance to unit sphere")
    plt.ylabel("Absolute error")
    plt.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(figure_path, dpi=200)
    print(f"\nSaved plot to {figure_path}")


if __name__ == "__main__":
    main()
