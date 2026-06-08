# Examples

The examples are organized by the scientific role they play in the package.
Files have not been moved yet, to avoid breaking notebook paths and existing
links.

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

- ` adapted_mesh_Gauss_Bonnet_benchmark.ipynb`

  Demonstrates a two-stage workflow: first construct a conforming
  indicator-refined reference mesh, then run a polynomial-degree convergence
  study on that adapted mesh.

## B. Singular And Near-Singular Integral Operators

These examples are experimental. They test square-squeezed product integration
(SSPI), Chebyshev-tail diagnostics, curvature-corrected singular models, and
hybrid SSPI-QBX ideas on analytic sphere problems.

- `laplace_single_layer_sphere_error.py`

  Tests singular, near-singular, and regular Laplace single-layer evaluations
  for constant density on the unit sphere.

- `laplace_single_layer_near_singular_convergence.py`

  Compares ordinary high-order quadrature with the SSPI prototype as the target
  approaches the surface.

- `laplace_single_layer_spherical_harmonic_error.py`

  Tests the Laplace single-layer prototype for a nonconstant degree-one
  spherical-harmonic density.

- `laplace_single_layer_moment_convergence.py`

  Varies the Chebyshev smooth degree and singular moment order to study
  convergence of the SSPI correction.

- `laplace_single_layer_diagnostics.py`

  Demonstrates base/enriched SSPI diagnostics, Chebyshev-tail indicators, and
  corrected-panel counts.

- `laplace_single_layer_hybrid_qbx.py`

  Tests the experimental hybrid SSPI-QBX evaluation path and prints QBX radius,
  order, convergence ratio, and diagnostic error estimates.

## C. PDE And Integral-Equation Prototypes

These examples are prototypes for future surface PDE and integral-equation
workflows. They are intentionally small and use manufactured or sphere-specific
solutions so that error can be measured.

- `harmonic_extension_integral_equation.py`

  Assembles a dense single-layer boundary integral equation for a harmonic
  extension problem on the unit sphere. This is a demonstration, not a reusable
  integral-equation solver API.

- `screened_laplace_beltrami_integral_operator.py`

  Applies the exact screened Green operator on the unit sphere to a
  manufactured solution of `(alpha - Delta_Gamma) u = f`.

- `screened_laplace_beltrami_parametrix.py`

  Studies the split of the screened Green kernel into a curvature-corrected
  singular parametrix and a smooth Chebyshev remainder.

- `screened_laplace_beltrami_sspi_parametrix.py`

  Connects the screened parametrix split to the SSPI quadrature path and
  compares against a known sphere solution.

- `pde/screened_laplace_beltrami_sphere_convergence.py`

  Benchmarks the screened Laplace-Beltrami prototype on the unit sphere using
  the manufactured solution `u(x,y,z)=z`. It separates geometry refinement,
  quadrature refinement, and smooth-remainder kernel approximation effects as
  far as the current prototype permits. This is not a general surface PDE
  solver.

## Notes

- The notebooks preserve the historical benchmark style and figures.
- The `.py` examples are easier to run in automated checks and are preferred
  for new prototype experiments.
- The singular and PDE examples should be treated as research implementations
  until the APIs and convergence behavior are broadened.
