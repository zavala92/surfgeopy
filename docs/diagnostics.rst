Diagnostics
===========

Integration Diagnostics
-----------------------

``integrate_with_diagnostics`` estimates numerical accuracy by comparing two
integration runs. The first run uses the user-provided ``IntegrationConfig``.
The second run enriches the interpolation and quadrature degree. The difference
between the two totals is reported as an a posteriori error estimate.

.. code-block:: python

   from surfgeopy import integrate_with_diagnostics

   report = integrate_with_diagnostics(
       surface,
       integrand,
       config,
       interpolation_degree_step=2,
       integration_degree_step=2,
       relative_tolerance=1.0e-8,
   )

   print(report.summary())

The returned ``DiagnosticIntegrationResult`` contains:

``total``
   The total integral from the enriched run.

``base_result``
   The original ``IntegrationResult``.

``enriched_result``
   The ``IntegrationResult`` computed with the enriched configuration.

``absolute_error_estimate``
   Absolute difference between enriched and base totals.

``relative_error_estimate``
   Absolute error estimate divided by the magnitude of the enriched total.

``local_absolute_errors``
   Per-face differences between enriched and base integral contributions.

``recommended_config``
   The next configuration to try if the requested tolerance was not reached.

``summary()``
   A compact text report suitable for notebooks, logs, and benchmark output.

For example, a report may look like:

.. code-block:: text

   Integral:              12.56637061435917
   Base integral:         12.56637061438142
   Estimated abs. error:  2.225e-11
   Estimated rel. error:  1.771e-12
   Max local error:       4.380e-13
   Quadrature points:     4704
   Interpolation degree:  8
   Integration degree:    16
   Recommendation:        accuracy target reached

The estimate is not a rigorous proof of the error. It is a practical numerical
diagnostic: if enriching the geometry and quadrature hardly changes the answer,
the computed integral is usually stable with respect to those parameters.

Adaptive Refinement
-------------------

``adaptive_integrate`` uses the local error indicators from
``integrate_with_diagnostics`` to refine only the faces that contribute most to
the estimated error.

.. code-block:: python

   from surfgeopy import adaptive_integrate

   adaptive = adaptive_integrate(
       surface,
       integrand,
       config,
       relative_tolerance=1.0e-8,
       max_iterations=4,
       marking_fraction=0.25,
   )

   print(adaptive.summary())

At each iteration, ``surfgeopy``:

1. Computes base and enriched integrals on the current mesh.
2. Forms local indicators from the per-face differences.
3. Marks the largest indicators according to ``marking_fraction``.
4. Refines the marked faces by triangular quadrisection.
5. Repeats until the tolerance is reached or ``max_iterations`` is exhausted.

The returned ``AdaptiveIntegrationResult`` contains:

``total``
   Final enriched integral value.

``absolute_error_estimate`` and ``relative_error_estimate``
   Final global diagnostic estimates.

``final_surface``
   The final ``LevelSetSurface`` with the adaptively refined reference mesh.

``history``
   Per-iteration records with face counts, marked faces, and error estimates.

``converged``
   Whether the requested diagnostic tolerance was reached.

The first implementation uses local triangular quadrisection. This may produce
hanging nodes in the reference mesh, but each face is integrated independently
and the refined faces still partition the original reference triangles.

Projection Diagnostics
----------------------

Closest-point projection is central to ``surfgeopy``. The package therefore
provides a diagnostic projection method in addition to the simple point-returning
method.

Simple Projection
-----------------

.. code-block:: python

   from surfgeopy import ImplicitSurface

   surface = ImplicitSurface(phi, grad_phi)
   projected = surface.project(point)

Projection With Diagnostics
---------------------------

.. code-block:: python

   info = surface.project_with_info(point)

   print(info.point)
   print(info.converged)
   print(info.iterations)
   print(info.residual)

The returned ``ProjectionResult`` contains:

``point``
   The projected point.

``converged``
   Whether the residual reached the configured tolerance.

``iterations``
   Number of projection iterations used.

``residual``
   Absolute value of ``phi(point)`` after projection.

Why This Matters
----------------

Projection failures usually indicate that the initial reference mesh is too far
from the implicit surface, the level-set gradient is near zero, or the maximum
iteration count is too small. Inspecting projection diagnostics helps identify
when integration error is caused by geometry projection rather than quadrature.
