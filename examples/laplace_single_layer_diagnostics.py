"""Chebyshev-tail diagnostics for experimental SSPI quadrature.

This example tests the base/enriched diagnostic path for singular and
near-singular Laplace single-layer evaluation on the unit sphere.  It prints
true errors, estimated local corrections, tail indicators, and panel counts.
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

    direction = mesh.vertices[0] / np.linalg.norm(mesh.vertices[0])
    radii = np.array([1.0, 1.001, 1.005, 1.01, 1.05])
    targets = radii[:, None] * direction[None, :]
    exact = exact_single_layer_unit_sphere(targets)

    config = SingularIntegrationConfig(
        interpolation_degree=6,
        regular_order=12,
        smooth_degree=6,
        moment_order=28,
        correction_order=14,
        near_threshold=1.5,
        singular_model="curvature",
    )
    operator = LaplaceSingleLayerOperator(surface, config)
    diagnostic = operator.evaluate_with_diagnostics(targets, target_tolerance=1.0e-8)
    error = np.abs(diagnostic.values - exact)

    print(
        "radius | value          | exact          | true error   | estimate     | tail sum     "
        "| max tail    | corrected | near | singular"
    )
    for index, radius in enumerate(radii):
        print(
            f"{radius:6.3f} | "
            f"{diagnostic.values[index]:.12e} | "
            f"{exact[index]:.12e} | "
            f"{error[index]:.3e} | "
            f"{diagnostic.absolute_error_estimates[index]:.3e} | "
            f"{diagnostic.tail_indicators[index]:.3e} | "
            f"{diagnostic.max_tail_indicators[index]:.3e} | "
            f"{diagnostic.corrected_panel_counts[index]:9d} | "
            f"{diagnostic.near_panel_counts[index]:4d} | "
            f"{diagnostic.singular_panel_counts[index]:8d}"
        )

    print("\nRecommended interpolation degree:", diagnostic.recommended_config.interpolation_degree)
    print("Recommended smooth degree:", diagnostic.recommended_config.smooth_degree)
    print("Recommended moment order:", diagnostic.recommended_config.moment_order)
    print("Recommended correction order:", diagnostic.recommended_config.correction_order)
    print("\nBase patches:", diagnostic.n_base_patches)
    print("Enriched patches built:", diagnostic.n_enriched_patches_built)
    print("Enriched patch fraction:", f"{diagnostic.enriched_patch_fraction:.3f}")
    print("Total corrected panel uses:", int(np.sum(diagnostic.corrected_panel_counts)))


if __name__ == "__main__":
    main()
