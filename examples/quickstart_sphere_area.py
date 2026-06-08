"""Quickstart: integrate 1 over the unit sphere."""

import numpy as np

from surfgeopy import IntegrationConfig, LevelSetSurface, integrate


def main() -> None:
    surface = LevelSetSurface.unit_sphere(mesh_refinement_level=1)
    config = IntegrationConfig(
        interpolation_degree=4,
        integration_degree=8,
        quadrature_rule="Gauss_Legendre",
    )

    result = integrate(surface, lambda _: 1.0, config)
    exact = 4.0 * np.pi

    print(f"area        : {result.total:.12f}")
    print(f"exact 4*pi  : {exact:.12f}")
    print(f"abs error   : {abs(result.total - exact):.3e}")
    print(f"quad points : {result.n_quadrature_points}")


if __name__ == "__main__":
    main()
