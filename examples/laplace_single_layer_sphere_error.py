"""Error check for prototype singular surface integration on the unit sphere.

For the Laplace single-layer kernel

    G(x, y) = 1 / (4*pi*|x-y|)

and constant density sigma = 1 on the unit sphere, the exact potential is

    u(x) = 1            for |x| <= 1,
    u(x) = 1 / |x|      for |x| > 1.

This gives a compact application where singular, near-singular, and regular
targets can all be checked against a known answer.
"""

from pathlib import Path

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

    surface_target = mesh.vertices[0] / np.linalg.norm(mesh.vertices[0])
    targets = np.array(
        [
            surface_target,
            1.01 * surface_target,
            1.10 * surface_target,
            np.array([0.0, 0.0, 2.0]),
        ]
    )

    configs = [
        (
            "prototype-low",
            SingularIntegrationConfig(
                interpolation_degree=4,
                regular_order=8,
                smooth_degree=4,
                moment_order=18,
                correction_order=10,
                near_threshold=1.5,
                singular_model="curvature",
            ),
        ),
        (
            "prototype-high",
            SingularIntegrationConfig(
                interpolation_degree=6,
                regular_order=12,
                smooth_degree=6,
                moment_order=28,
                correction_order=14,
                near_threshold=1.5,
                singular_model="curvature",
            ),
        ),
    ]

    exact = exact_single_layer_unit_sphere(targets)

    for label, config in configs:
        operator = LaplaceSingleLayerOperator(surface, config)
        result = operator.evaluate(targets)
        error = np.abs(result.values - exact)

        print(f"\n{label}")
        print(
            "target | |x|      | numerical      | exact          | abs error      "
            "| near | singular"
        )
        for index, target in enumerate(targets):
            print(
                f"{index:6d} | "
                f"{np.linalg.norm(target):.6f} | "
                f"{result.values[index]:.12e} | "
                f"{exact[index]:.12e} | "
                f"{error[index]:.3e} | "
                f"{result.near_panel_counts[index]:4d} | "
                f"{result.singular_panel_counts[index]:8d}"
            )


if __name__ == "__main__":
    main()
