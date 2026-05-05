# Surfgeopy Roadmap

The goal is to turn `surfgeopy` from a research prototype into a reliable, benchmarked library for high-order surface integration on embedded surfaces.

## 1. Baseline Reliability

- Keep `pytest` green from a fresh checkout.
- Remove generated operating-system files from source distributions.
- Align Python version metadata across packaging files.
- Add continuous integration for supported Python versions.
- Add deterministic tests for analytic integrals on canonical surfaces.

## 2. Numerical Kernel

- Separate geometry projection, geometry interpolation, quadrature generation, and integration accumulation.
- Add clear interfaces for implicit surfaces, explicit meshes, and future CAD-style surface patches.
- Support vectorized scalar and vector-valued integrands.
- Add automatic quadrature degree and refinement selection.
- Add robust handling for near-degenerate triangles and high-curvature patches.

## 3. Benchmarks

- Build reproducible benchmark cases for spheres, ellipsoids, tori, genus-two surfaces, bioconcave surfaces, and Dziuk surfaces.
- Track error, convergence rate, runtime, memory use, and failure modes.
- Compare against reference methods for implicit-surface and high-order embedded-boundary quadrature.
- Publish benchmark scripts and result tables with fixed random seeds and mesh versions.

## 4. Public API

- Introduce high-level objects such as `SurfaceMesh`, `LevelSetSurface`, and `QuadratureRule`.
- Provide a simple `integrate(surface, integrand, ...)` entry point.
- Keep lower-level kernels available for research experiments.
- Add type hints and API documentation for all public functions.

## 5. Documentation And Releases

- Add installation instructions for users and contributors.
- Convert notebooks into tested examples where possible.
- Publish API docs and benchmark documentation.
- Produce versioned releases with changelogs and source distributions.

