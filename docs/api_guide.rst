API Guide
=========

Modern API
----------

The main user-facing classes are:

``SurfaceMesh``
   Stores reference mesh vertices and faces. Use ``SurfaceMesh.icosphere`` for
   a built-in sphere mesh, ``SurfaceMesh.from_mat`` for MATLAB mesh files, or
   ``SurfaceMesh.from_gmsh`` / ``SurfaceMesh.to_gmsh`` for ASCII Gmsh surface
   meshes. High-order triangular Gmsh elements are linearized by default; pass
   ``preserve_order=True`` when loading to keep their full node connectivity.

``LevelSetSurface``
   Bundles a ``SurfaceMesh`` with ``phi`` and ``grad_phi``. For examples and
   first experiments, ``LevelSetSurface.unit_sphere`` and
   ``LevelSetSurface.sphere`` create a ready-to-use spherical surface.

``IntegrationConfig``
   Stores interpolation degree, refinement level, integration degree,
   quadrature rule, and the polynomial ``lp`` degree.

``IntegrationResult``
   Stores per-face values, quadrature points, weights, offsets, and the total.

``CompiledIntegrationScheme``
   Stores reusable quadrature points, weights, offsets, and configuration for
   repeated integrations on the same surface.

``DiagnosticIntegrationResult``
   Stores two integration runs, an a posteriori error estimate, local per-face
   differences, and a text summary.

``SurfaceGeometryResult``
   Stores quadrature-point samples of the interpolated surface geometry:
   tangents, normal, metric tensor, area density, second fundamental form, mean
   curvature, and Gaussian curvature.

``IndicatorRefinementResult``
   Stores the final surface and per-step history from indicator-based reference
   mesh refinement, including indicator statistics, marked-face fractions,
   mesh growth counts, and the terminal stop reason.

``integrate``
   Runs the high-order implicit-surface integration workflow.

``compile_integration``
   Builds reusable quadrature once. Use this for parameter sweeps or many
   integrands over a fixed surface and configuration.

``integrate_with_diagnostics``
   Runs a base and enriched integration configuration to estimate accuracy.

``refine_by_indicator``
   Refines the reference mesh using an indicator evaluated at affine face
   centers. This is useful when the adapted mesh should be generated before a
   separate polynomial degree study.

``surface_geometry``
   Evaluates Minterpy spectral derivatives of the high-order surface map and
   returns differential geometry quantities at the quadrature points.

Experimental Singular And PDE Prototype APIs
--------------------------------------------

The following names are available from ``surfgeopy`` but should be treated as
research implementations. They are useful for reproducing the included
singular-kernel and screened Laplace-Beltrami experiments, but they are not yet
a stable general integral-equation API.

``SingularIntegrationConfig``
   Configuration for experimental square-squeezed product integration (SSPI)
   of Laplace single-layer potentials.

``LaplaceSingleLayerOperator``
   Reusable prototype evaluator for the Laplace single-layer potential. It
   supports an SSPI path and an experimental hybrid SSPI-QBX path.

``SingularIntegralResult``
   Potential values and near/singular panel counts. Hybrid QBX runs also report
   QBX flags, radii, orders, convergence ratios, and diagnostic estimates.

``SingularDiagnosticResult``
   Base/enriched diagnostic data for SSPI, including Chebyshev-tail indicators
   and corrected-panel counts.

``ScreenedParametrixConfig``
   Configuration for the screened Laplace-Beltrami parametrix prototype.

``ScreenedLaplaceBeltramiParametrixOperator``
   Prototype evaluator for a split screened Green kernel: local singular
   parametrix plus a user-supplied smooth remainder.

``ScreenedParametrixResult``
   Potential values, near/singular panel counts, and accumulated tail
   indicators for the screened parametrix prototype.

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

Repeated Integrations
---------------------

.. code-block:: python

   scheme = compile_integration(surface, config)

   area = scheme.integrate(lambda _: 1.0)
   moment = scheme.integrate(lambda x: x[2] ** 2)

   vectorized_moment = scheme.integrate(
       lambda points: points[:, 2] ** 2,
       vectorized=True,
   )

``compile_integration`` avoids rebuilding projected high-order geometry when
only the scalar field changes. With ``vectorized=True``, the integrand receives
all quadrature points at once as an array with shape ``(n_points, 3)`` and must
return one scalar per row.

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

Indicator Refinement Example
----------------------------

.. code-block:: python

   adapted = refine_by_indicator(
       surface,
       integrand,
       max_iterations=6,
       threshold_fraction=0.25,
   )

   adapted_surface = adapted.final_surface

   errors = []
   for degree in range(2, 22):
       config = IntegrationConfig(
           interpolation_degree=degree,
           integration_degree=15,
           quadrature_rule="Pull_back_Gauss",
       )
       result = integrate(adapted_surface, integrand, config)
       errors.append(abs(result.total - exact_value) / abs(exact_value))

``refine_by_indicator`` follows the host-grid adaptation pattern: it evaluates
``abs(indicator(center))`` on each linear reference triangle, marks faces above
``threshold_fraction * max_indicator``, and subdivides the marked faces. By
default, adjacent faces with split edges are also green-refined so the adapted
reference mesh remains conforming. The high-order curved interpolation is then
constructed on the adapted reference mesh during the subsequent integration
run.

The returned ``IndicatorRefinementResult`` keeps the full per-iteration
diagnostic trace in ``history``. Convenience arrays such as ``face_counts``,
``vertex_counts``, ``max_indicators``, ``mean_indicators``, and
``marked_fractions`` are useful for logging or plotting adaptation behavior.
``stop_reason`` reports whether the run stopped because the indicator vanished,
no faces were above threshold, the mesh was empty, or ``max_iterations`` was
reached.

Surface Geometry Example
------------------------

.. code-block:: python

   geometry = surface_geometry(surface, config)

   print(geometry.points)
   print(geometry.normal)
   print(geometry.metric_tensor)
   print(geometry.area_density)
   print(geometry.gaussian_curvature)

``surface_geometry`` uses the same cubical reparametrization and Minterpy
interpolation backend as ``integrate``. The interpolant is built on Minterpy's
Chebyshev-Lobatto grid, and the resulting Newton polynomial is differentiated
with respect to the two square coordinates. First derivatives give the surface
tangents and metric tensor; second derivatives give the second fundamental form
and curvatures. Curvature computation requires ``interpolation_degree >= 2``.

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
