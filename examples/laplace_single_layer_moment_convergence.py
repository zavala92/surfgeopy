"""Moment-order convergence for the prototype SSPI quadrature.

This keeps the surface mesh, targets, and geometry degree fixed, then increases
the Chebyshev smooth degree and singular moment quadrature order. The exact
Laplace single-layer potential is known on and outside the unit sphere.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from surfgeopy import (
    LaplaceSingleLayerOperator,
    LevelSetSurface,
    SingularIntegrationConfig,
    SurfaceMesh,
)


def phi(x: np.ndarray) -> float:
    return x[0] ** 2 + x[1] ** 2 + x[2] ** 2 - 1.0


def dphi(x: np.ndarray) -> np.ndarray:
    return np.array([2.0 * x[0], 2.0 * x[1], 2.0 * x[2]])


def exact_single_layer_unit_sphere(targets: np.ndarray) -> np.ndarray:
    radii = np.linalg.norm(targets, axis=1)
    return np.where(radii <= 1.0, 1.0, 1.0 / radii)


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    mesh = SurfaceMesh.from_mat(str(root / "meshes" / "SphereMesh_N=124_r=1.mat"))
    surface = LevelSetSurface(mesh, phi, dphi)

    direction = mesh.vertices[0] / np.linalg.norm(mesh.vertices[0])
    radii = np.array([1.0, 1.001, 1.01])
    targets = radii[:, None] * direction[None, :]
    exact = exact_single_layer_unit_sphere(targets)

    configurations = [
        (3, 14),
        (4, 18),
        (5, 24),
        (6, 32),
        (8, 44),
    ]

    errors = []
    print("smooth degree | moment order | surface error | r=1.001 error | r=1.01 error")
    for smooth_degree, moment_order in configurations:
        config = SingularIntegrationConfig(
            interpolation_degree=6,
            regular_order=12,
            smooth_degree=smooth_degree,
            moment_order=moment_order,
            correction_order=max(12, smooth_degree + 6),
            near_threshold=1.5,
            singular_model="curvature",
        )
        operator = LaplaceSingleLayerOperator(surface, config)
        result = operator.evaluate(targets)
        error = np.abs(result.values - exact)
        errors.append(error)
        print(
            f"{smooth_degree:13d} | "
            f"{moment_order:12d} | "
            f"{error[0]:13.6e} | "
            f"{error[1]:13.6e} | "
            f"{error[2]:12.6e}"
        )

    errors = np.asarray(errors)
    labels = ["surface", "r=1.001", "r=1.01"]

    figure_path = root / "images" / "laplace_single_layer_moment_convergence.png"
    figure_path.parent.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(7.0, 4.8))
    for target_index, label in enumerate(labels):
        plt.semilogy(
            [item[0] for item in configurations],
            errors[:, target_index],
            "-o",
            label=label,
        )
    plt.xlabel("Chebyshev smooth degree")
    plt.ylabel("Absolute error")
    plt.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(figure_path, dpi=200)
    print(f"\nSaved plot to {figure_path}")


if __name__ == "__main__":
    main()
