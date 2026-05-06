API Guide
=========

Modern API
----------

The main user-facing classes are:

``SurfaceMesh``
   Stores reference mesh vertices and faces. Use ``SurfaceMesh.from_mat`` for
   the MATLAB mesh files distributed with the examples.

``LevelSetSurface``
   Bundles a ``SurfaceMesh`` with ``phi`` and ``grad_phi``.

``IntegrationConfig``
   Stores interpolation degree, refinement level, integration degree,
   quadrature rule, and the polynomial ``lp`` degree.

``IntegrationResult``
   Stores per-face values, quadrature points, weights, offsets, and the total.

``DiagnosticIntegrationResult``
   Stores two integration runs, an a posteriori error estimate, local per-face
   differences, and a text summary.

``integrate``
   Runs the high-order implicit-surface integration workflow.

``integrate_with_diagnostics``
   Runs a base and enriched integration configuration to estimate accuracy.

``adaptive_integrate``
   Repeats diagnostic integration, refines faces with the largest local error
   indicators, and returns convergence history.

Example
-------

.. code-block:: python

   config = IntegrationConfig(
       interpolation_degree=8,
       refinement_level=1,
       integration_degree=14,
       quadrature_rule="Pull_back_Gauss",
   )

   result = integrate(surface, integrand, config)

Diagnostics Example
-------------------

.. code-block:: python

   report = integrate_with_diagnostics(
       surface,
       integrand,
       config,
       interpolation_degree_step=2,
       integration_degree_step=2,
       relative_tolerance=1.0e-8,
   )

   print(report.total)
   print(report.absolute_error_estimate)
   print(report.relative_error_estimate)
   print(report.summary())

``report.total`` is the value from the enriched run. The error estimate is the
difference between the enriched and base totals. If the requested tolerance is
not reached, ``report.recommended_config`` proposes the next enriched
configuration to try.

Adaptive Example
----------------

.. code-block:: python

   adaptive = adaptive_integrate(
       surface,
       integrand,
       config,
       relative_tolerance=1.0e-8,
       max_iterations=4,
       marking_fraction=0.25,
   )

   print(adaptive.total)
   print(adaptive.relative_error_estimate)
   print(adaptive.summary())

``adaptive_integrate`` marks the faces with the largest local error indicators
and refines them by triangular quadrisection. The result stores the final
surface mesh and an iteration history with the number of faces, marked faces,
and error estimate at each step.

Choosing Degrees
----------------

``interpolation_degree``
   Controls the polynomial degree used to approximate the curved geometry.
   Higher values can improve accuracy but increase cost.

``integration_degree``
   Controls the quadrature rule. Higher values improve integration accuracy for
   the selected geometry approximation.

``refinement_level``
   Applies triangular quadrisection before constructing curved patches. This is
   useful when the initial mesh is coarse or the surface curvature is high.

``lp_degree``
   Controls the polynomial multi-index degree in the interpolation backend.
   ``float("inf")`` selects tensor-product-style degree behavior.

Validation
----------

``IntegrationConfig`` validates that interpolation and integration degrees are
positive and that refinement is non-negative.
