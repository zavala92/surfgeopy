"""Experimental hybrid SSPI-QBX check for a Laplace layer potential.

This example compares the SSPI path with the prototype QBX near-surface
evaluation path for the Laplace single-layer potential on the unit sphere.  It
is intended as a diagnostic experiment, not as a general QBX solver.
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
    radii = np.array([1.0, 1.001, 1.01, 1.10, 2.00])
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
    hybrid_qbx_config = SingularIntegrationConfig(
        interpolation_degree=6,
        regular_order=12,
        smooth_degree=6,
        moment_order=28,
        correction_order=14,
        near_threshold=1.5,
        singular_model="curvature",
        evaluation_strategy="hybrid_qbx",
        qbx_order=12,
        qbx_quadrature_order=36,
        qbx_radius_factor=0.5,
        qbx_radius_padding=0.25,
    )

    sspi = LaplaceSingleLayerOperator(surface, sspi_config).evaluate(targets)
    hybrid_qbx = LaplaceSingleLayerOperator(surface, hybrid_qbx_config).evaluate(targets)

    print(
        "radius | exact          | SSPI error    | QBX error     | QBX est.     | "
        "rho    | qbx r. | near | singular | QBX"
    )
    for index, radius in enumerate(radii):
        sspi_error = abs(sspi.values[index] - exact[index])
        qbx_error = abs(hybrid_qbx.values[index] - exact[index])
        print(
            f"{radius:6.3f} | "
            f"{exact[index]:.12e} | "
            f"{sspi_error:.3e} | "
            f"{qbx_error:.3e} | "
            f"{hybrid_qbx.qbx_error_estimates[index]:.3e} | "
            f"{hybrid_qbx.qbx_convergence_ratios[index]:.3f} | "
            f"{hybrid_qbx.qbx_radii[index]:.3e} | "
            f"{hybrid_qbx.near_panel_counts[index]:4d} | "
            f"{hybrid_qbx.singular_panel_counts[index]:8d} | "
            f"{bool(hybrid_qbx.qbx_target_flags[index])}"
        )


if __name__ == "__main__":
    main()
