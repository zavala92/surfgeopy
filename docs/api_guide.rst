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

``integrate``
   Runs the high-order implicit-surface integration workflow.

Example
-------

.. code-block:: python

   config = IntegrationConfig(
       interpolation_degree=8,
       lp_degree=float("inf"),
       refinement_level=1,
       integration_degree=14,
       quadrature_rule="Pull_back_Gauss",
   )

   result = integrate(surface, integrand, config)

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

