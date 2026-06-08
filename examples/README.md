# Examples

The examples are organized by the scientific role they play in the package.
They currently focus on smooth high-order surface quadrature.

## A. Core Geometry And Quadrature

These examples exercise the mature part of `surfgeopy`: high-order
square-squeezed quadrature on implicitly represented curved surface
triangulations.

### First Run

- `quickstart_sphere_area.py`

  Minimal script using `LevelSetSurface.unit_sphere(...)` and the built-in
  icosphere reference mesh. This is the best first file to run after installing
  the package.

### Smooth Surface Integration

- `Test_integration_on_the_whole_sphere.ipynb`

  Tests high-order convergence for sphere area and a spherical-harmonic
  integrand with a known zero integral.

- `Test_integration_on_torus.ipynb`

  Tests high-order surface-area integration on a torus represented by a
  level-set function.

### Curvature And Gauss-Bonnet Benchmarks

- `Gauss_Bonnet_theorem_bench.ipynb`

  Tests high-order integration of Gaussian-curvature expressions on surfaces
  with known Euler characteristic.

- `adapted_mesh_Gauss_Bonnet_benchmark.ipynb`

  Demonstrates a two-stage workflow: first construct a conforming
  indicator-refined reference mesh, then run a polynomial-degree convergence
  study on that adapted mesh.

## Notes

- The notebooks preserve the historical benchmark style and figures.
- The `.py` examples are easier to run in automated checks and are preferred
  for new examples.
