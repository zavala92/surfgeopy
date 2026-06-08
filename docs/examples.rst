Examples
========

The examples are grouped by their role in the numerical story. The regular
surface-quadrature examples are the mature part of the package. Singular,
near-singular, and PDE-oriented examples are research prototypes unless stated
otherwise.

Core Geometry And Quadrature
----------------------------

These examples test high-order integration on curved triangulated implicit
surfaces using square-squeezed geometry interpolation and high-order
quadrature.

First script to run:

``examples/quickstart_sphere_area.py``
   Minimal unit-sphere area example using ``LevelSetSurface.unit_sphere`` and
   the built-in icosphere reference mesh.

.. toctree::
   :maxdepth: 1

   examples/sphere
   examples/torus
   examples/gauss_bonnet

Notebook examples in the repository:

``examples/Test_integration_on_the_whole_sphere.ipynb``
   Sphere-area and spherical-harmonic benchmarks.

``examples/Test_integration_on_torus.ipynb``
   Torus-area convergence benchmark.

``examples/Gauss_Bonnet_theorem_bench.ipynb``
   Gauss-Bonnet curvature-integral benchmarks.

``examples/ adapted_mesh_Gauss_Bonnet_benchmark.ipynb``
   Indicator-refined Gauss-Bonnet benchmark on a biconcave surface.

Singular And Near-Singular Operators
------------------------------------

These examples test experimental square-squeezed product integration (SSPI),
Chebyshev-tail diagnostics, curvature-corrected singular models, and hybrid
SSPI-QBX evaluation.

``examples/laplace_single_layer_sphere_error.py``
   Laplace single-layer error check for singular, near-singular, and regular
   targets on the unit sphere.

``examples/laplace_single_layer_near_singular_convergence.py``
   Comparison between ordinary surface quadrature and SSPI as targets approach
   the surface.

``examples/laplace_single_layer_spherical_harmonic_error.py``
   Nonconstant-density Laplace single-layer test with a degree-one spherical
   harmonic.

``examples/laplace_single_layer_moment_convergence.py``
   Moment-order and smooth-degree convergence for SSPI.

``examples/laplace_single_layer_diagnostics.py``
   Base/enriched diagnostics, Chebyshev-tail indicators, and corrected-panel
   counts.

``examples/laplace_single_layer_hybrid_qbx.py``
   Experimental hybrid SSPI-QBX near-surface evaluation.

PDE And Integral-Equation Prototypes
------------------------------------

These examples are early demonstrations for future surface integral-equation
workflows. They use sphere-specific manufactured solutions so the error can be
measured.

``examples/harmonic_extension_integral_equation.py``
   Dense single-layer boundary-integral demonstration for harmonic extension
   on the unit sphere.

``examples/screened_laplace_beltrami_integral_operator.py``
   Exact spherical screened Green operator applied to
   ``(alpha - Delta_Gamma) u = f``.

``examples/screened_laplace_beltrami_parametrix.py``
   Curvature-corrected screened singular parametrix plus smooth Chebyshev
   remainder.

``examples/screened_laplace_beltrami_sspi_parametrix.py``
   SSPI evaluation of the screened parametrix split on the unit sphere.

``examples/pde/screened_laplace_beltrami_sphere_convergence.py``
   Screened Laplace-Beltrami benchmark with manufactured solution
   ``u(x,y,z)=z``. The benchmark reports geometry-refinement errors,
   quadrature-refinement errors, and a kernel-level smooth-remainder study.
   It is a first step toward surface PDE benchmarks, not a general solver.
