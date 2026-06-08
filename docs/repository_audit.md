# Repository Audit

This audit records the current implementation state of `surfgeopy` relative to
the intended package story:

> high-order square-squeezed quadrature for regular, singular, and
> near-singular integral operators on implicitly represented curved surface
> triangulations.

The purpose is to distinguish implemented functionality from experimental
modules, prototypes, and future work.

## Implemented Core

### Implicit Surface Geometry

Implemented in:

- `surfgeopy/surface.py`
- `surfgeopy/api.py`
- `surfgeopy/surf_integration.py`
- `surfgeopy/utils.py`

Current functionality:

- `ImplicitSurface` stores a level-set function and gradient.
- `project_with_info` performs closest-point projection and reports convergence
  diagnostics.
- `LevelSetSurface` combines a reference mesh with `phi` and `grad_phi`.
- `surface_geometry` differentiates the high-order Minterpy surface map and
  returns tangents, normals, metric tensors, area density, second fundamental
  form, mean curvature, and Gaussian curvature.

Status: implemented and tested on sphere cases.

### Triangulated Curved Surfaces

Implemented in:

- `surfgeopy/api.py`
- `surfgeopy/surface.py`
- `surfgeopy/surf_integration.py`
- `surfgeopy/remesh.py`

Current functionality:

- `SurfaceMesh` stores vertices and triangular or polygonal faces.
- `.mat` meshes are loaded through `SurfaceMesh.from_mat`.
- Each face is decomposed into triangles when needed.
- Curved patches are built by projecting interpolation nodes from each affine
  reference triangle onto the implicit surface.
- `subdivide` and `subdivide_conforming` provide uniform and local conforming
  reference-mesh refinement.

Status: implemented. The conforming indicator refinement is useful for degree
studies but is still a simple red/green refinement strategy, not a complete
adaptive error-control framework.

### Square-Squeezed Parametrizations

Implemented in:

- `surfgeopy/utils.py`
- `surfgeopy/reference_quadrature.py`
- `surfgeopy/surf_integration.py`
- `surfgeopy/singular_integrals.py`

Current functionality:

- `pushforward` maps square interpolation nodes to the reference simplex.
- `pullback` maps simplex quadrature points back to the square.
- The non-Duffy branch implements the square-squeezing map used by the package.
- Simplex rules are evaluated through the square-squeezed geometry map by
  applying the appropriate weight scale.

Status: implemented and central to the package.

### Chebyshev-Lobatto Interpolation

Implemented through:

- `surfgeopy/surf_integration.py`
- `surfgeopy/singular_integrals.py`
- Minterpy's `MultiIndexSet`, `Grid`, `NewtonPolynomial`, and `dds`

Current functionality:

- Minterpy provides the interpolation grid and Newton polynomial backend.
- The surface map is interpolated on the square.
- First and second polynomial derivatives are used for quadrature weights,
  metric tensors, curvature, and curvature-corrected singular models.

Status: implemented. The package relies on Minterpy for the interpolation
backend.

### Smooth Surface Quadrature

Implemented in:

- `surfgeopy/api.py`
- `surfgeopy/surf_integration.py`
- `surfgeopy/reference_quadrature.py`
- `surfgeopy/quadrature_points.py`
- `surfgeopy/quadrature_points_gl.py`

Current functionality:

- `integrate` is the modern high-level API.
- `integration` is the legacy compatibility function.
- `IntegrationConfig` controls interpolation degree, quadrature degree,
  refinement level, `lp_degree`, and quadrature rule.
- `IntegrationResult` stores per-face values, quadrature points, weights,
  offsets, and total integral.
- `integrate_with_diagnostics` compares base and enriched configurations.
- Available quadrature families include `ModePy_VioreanuRokhlin`,
  `ModePy_XiaoGimbutas`, `ModePy_GrundmannMoeller`, `Pull_back_Gauss`, and
  `Gauss_Legendre`.

Status: implemented and tested.

Examples:

- `examples/Test_integration_on_the_whole_sphere.ipynb`
- `examples/Test_integration_on_torus.ipynb`
- `examples/Gauss_Bonnet_theorem_bench.ipynb`
- `examples/ adapted_mesh_Gauss_Bonnet_benchmark.ipynb`
- `benchmarks/benchmark_sphere_area.py`
- `benchmarks/run_benchmarks.py`

## Experimental Singular And Near-Singular Integration

### Singular Product Integration And SSPI

Implemented in:

- `surfgeopy/singular_integrals.py`

Current functionality:

- `SingularIntegrationConfig`
- `LaplaceSingleLayerOperator`
- `laplace_single_layer_potential`
- square-based Chebyshev interpolation of the smooth factor,
- Duffy-style moment integration for singular Chebyshev moments,
- regular quadrature away from near panels.

Status: experimental research implementation. It is validated on sphere
single-layer examples, but it is not yet a general boundary integral library.

Examples:

- `examples/laplace_single_layer_sphere_error.py`
- `examples/laplace_single_layer_moment_convergence.py`
- `examples/laplace_single_layer_spherical_harmonic_error.py`

### Near-Singular Evaluation

Implemented in:

- `surfgeopy/singular_integrals.py`

Current functionality:

- nearest-patch search by square-parameter optimization,
- near and singular panel classification,
- product integration on near panels,
- regular quadrature on far panels.

Status: experimental. The current target classification is geometric and local;
automatic tolerance selection remains future work.

Examples:

- `examples/laplace_single_layer_near_singular_convergence.py`
- `examples/laplace_single_layer_diagnostics.py`

### Chebyshev-Tail Diagnostics

Implemented in:

- `surfgeopy/singular_integrals.py`

Current functionality:

- `SingularDiagnosticResult`
- `LaplaceSingleLayerOperator.evaluate_with_diagnostics`
- base/enriched comparison on near and singular panels,
- accumulated and maximum Chebyshev-tail indicators,
- counts of corrected panels and lazily built enriched patches.

Status: experimental diagnostic tool. The indicators are practical numerical
signals, not rigorous error bounds.

Example:

- `examples/laplace_single_layer_diagnostics.py`

### Curvature-Corrected Singular Models

Implemented in:

- `surfgeopy/singular_integrals.py`

Current functionality:

- `singular_model="metric"` uses the local metric distance.
- `singular_model="curvature"` uses first and second derivatives of the
  Minterpy surface map to form a quadratic/quartic local distance model.
- The screened parametrix prototype also includes curvature-corrected chord
  distances using a Gaussian-curvature expansion.

Status: experimental. The implementation is designed for testing and
comparison, not yet for broad robustness claims.

Examples:

- `examples/laplace_single_layer_sphere_error.py`
- `examples/laplace_single_layer_near_singular_convergence.py`
- `examples/screened_laplace_beltrami_parametrix.py`
- `examples/screened_laplace_beltrami_sspi_parametrix.py`

### QBX And Hybrid SSPI-QBX

Implemented in:

- `surfgeopy/singular_integrals.py`

Current functionality:

- `evaluation_strategy="hybrid_qbx"` for `LaplaceSingleLayerOperator`.
- automatic local center and radius selection based on closest-patch data,
- Legendre expansion for the Laplace kernel around a QBX center,
- diagnostic order comparison, radius, order, and convergence-ratio reporting.

Status: experimental prototype. It is currently tied to the Laplace
single-layer kernel and small validation cases.

Example:

- `examples/laplace_single_layer_hybrid_qbx.py`

### Layer-Potential Kernels

Implemented in:

- `surfgeopy/singular_integrals.py`
- `examples/harmonic_extension_integral_equation.py`

Current functionality:

- Laplace single-layer potential evaluation with constant and nonconstant
  densities.
- A small direct Nyström-style harmonic extension example on the unit sphere.

Status: experimental. No general integral-equation assembly API is currently
provided.

## PDE And Integral-Equation Prototypes

### Harmonic Extension

Implemented as an example:

- `examples/harmonic_extension_integral_equation.py`

Current functionality:

- assembles a dense single-layer matrix on quadrature points,
- uses a simple disk self-panel correction,
- evaluates a sphere Dirichlet test with known density and exact
  interior/exterior harmonic extension.

Status: prototype example. It is not a reusable solver API.

### Screened Laplace-Beltrami Prototypes

Implemented in:

- `surfgeopy/singular_integrals.py`
- `examples/screened_laplace_beltrami_integral_operator.py`
- `examples/screened_laplace_beltrami_parametrix.py`
- `examples/screened_laplace_beltrami_sspi_parametrix.py`
- `examples/pde/screened_laplace_beltrami_sphere_convergence.py`

Current functionality:

- exact spectral Green operator on the unit sphere for manufactured solutions,
- local screened parametrix split using `K0(sqrt(alpha) r)/(2*pi)`,
- Chebyshev approximation of the smooth remainder on the sphere,
- SSPI evaluation of the singular screened parametrix plus supplied smooth
  remainder.
- a focused sphere convergence benchmark for the manufactured solution
  `u(x,y,z)=z`, separating geometry refinement, quadrature refinement, and
  kernel-level smooth-remainder approximation effects.

Status: prototype. These examples show a route toward
`(alpha - Delta_Gamma) u = f` integral-operator experiments, but they do not
implement a general screened Laplace-Beltrami solver on arbitrary implicit
surfaces.

## Tests And Validation Coverage

Current tests in `tests/test_surfgeopy.py` cover:

- square/simplex pullback and pushforward maps,
- implicit projection diagnostics,
- mesh loading,
- smooth sphere area integration,
- integration diagnostics,
- conforming indicator refinement,
- surface geometry and curvature on a sphere,
- ModePy quadrature rule selection,
- Laplace single-layer SSPI on the unit sphere,
- near-singular off-surface Laplace target evaluation,
- hybrid SSPI-QBX diagnostic data,
- screened Laplace-Beltrami parametrix prototype on the sphere.
- screened Laplace-Beltrami sphere convergence benchmark smoke test.

Benchmark infrastructure in `benchmarks/` covers smooth area/convergence cases.

## Missing Or Future Work

Missing as package-level functionality:

- a general surface PDE solver API,
- general integral-equation assembly with boundary conditions and solver
  choices,
- arbitrary-surface screened Green remainders,
- fast summation or acceleration for many-source/many-target evaluation,
- rigorous a posteriori error estimators for singular operators,
- automatic quadrature/order/radius selection with convergence guarantees,
- broad robustness handling for near-degenerate triangles and difficult
  projections,
- tested vector-valued kernels and double-layer operators,
- a stable public API for all experimental singular and PDE features.

Recommended interpretation:

- Smooth high-order quadrature is the current core package capability.
- Singular and near-singular operators are experimental research modules.
- Harmonic extension and screened Laplace-Beltrami files are prototypes used to
  test the numerical direction.
