Quickstart
==========

This page shows the recommended modern API for integrating a scalar function
over an implicit surface.

The basic workflow is:

1. Create or load a triangulated reference mesh.
2. Define the implicit surface through ``phi`` and ``grad_phi``.
3. Choose an ``IntegrationConfig``.
4. Call ``integrate`` and read ``result.total``.

Unit Sphere Area
----------------

The unit sphere is represented as the zero level set

.. math::

   \phi(x, y, z) = x^2 + y^2 + z^2 - 1.

Its exact area is :math:`4\pi`.

From a terminal, the same smoke test is available as:

.. code-block:: bash

   surfgeopy demo

.. code-block:: python

   import numpy as np

   from surfgeopy import IntegrationConfig, LevelSetSurface, integrate


   surface = LevelSetSurface.unit_sphere(mesh_refinement_level=1)
   config = IntegrationConfig(
       interpolation_degree=4,
       integration_degree=8,
       quadrature_rule="Gauss_Legendre",
   )

   result = integrate(surface, lambda _: 1.0, config)
   print(result.total)
   print(4.0 * np.pi)
   print(result.n_quadrature_points)

For custom surfaces, build a ``SurfaceMesh`` and pass your own level-set
functions:

.. code-block:: python

   import numpy as np

   from surfgeopy import LevelSetSurface, SurfaceMesh


   def phi(x: np.ndarray) -> float:
       return x[0] ** 2 + x[1] ** 2 + x[2] ** 2 - 1.0


   def grad_phi(x: np.ndarray) -> np.ndarray:
       return np.array([2.0 * x[0], 2.0 * x[1], 2.0 * x[2]])


   mesh = SurfaceMesh.icosphere(refinement_level=1)
   surface = LevelSetSurface(mesh, phi, grad_phi)

Result Object
-------------

``integrate`` returns an ``IntegrationResult`` with:

``values``
   One integral contribution per input mesh face.

``points``
   Physical quadrature points on the curved surface approximation.

``weights``
   Physical quadrature weights.

``offsets``
   Face-to-quadrature-point offsets. Points for face ``i`` are in
   ``points[offsets[i]:offsets[i + 1]]``.

``total``
   The total integral over the surface.

Faster Repeated Integrations
----------------------------

If the surface and configuration stay fixed, build the quadrature once and
reuse it for many scalar fields:

.. code-block:: python

   from surfgeopy import compile_integration

   scheme = compile_integration(surface, config)
   area = scheme.integrate(lambda _: 1.0)
   z_second_moment = scheme.integrate(lambda x: x[2] ** 2)

For NumPy-style functions, evaluate all quadrature points in one call:

.. code-block:: python

   z_second_moment = scheme.integrate(
       lambda points: points[:, 2] ** 2,
       vectorized=True,
   )

Estimate Accuracy
-----------------

For accuracy-sensitive runs, use ``integrate_with_diagnostics`` to compare the
selected configuration with an enriched configuration:

.. code-block:: python

   from surfgeopy import integrate_with_diagnostics

   report = integrate_with_diagnostics(surface, lambda _: 1.0, config)
   print(report.summary())

The reported value ``report.total`` is the enriched integral. The fields
``absolute_error_estimate`` and ``relative_error_estimate`` measure the
difference between the base and enriched runs.

Surface Geometry
----------------

The high-order surface map can also be differentiated through Minterpy. The
geometry interpolant is built on the Chebyshev-Lobatto grid and differentiated
as a polynomial, which gives tangents, unit normals, metric tensors, area
densities, second fundamental forms, and curvatures at the quadrature points:

.. code-block:: python

   from surfgeopy import surface_geometry

   geometry = surface_geometry(surface, config)
   print(geometry.normal)
   print(geometry.gaussian_curvature)
   print(geometry.mean_curvature)

For the unit sphere, the Gaussian curvature is close to one and the absolute
mean curvature is close to one. Curvature quantities require
``interpolation_degree >= 2`` because second derivatives of the surface
interpolant are used.

Indicator-Based Mesh Adaptation
-------------------------------

For workflows where the mesh is adapted first and the polynomial degree is
varied afterwards, use ``refine_by_indicator``:

.. code-block:: python

   from surfgeopy import refine_by_indicator

   adapted = refine_by_indicator(
       surface,
       integrand,
       max_iterations=6,
       threshold_fraction=0.25,
   )

   adapted_surface = adapted.final_surface

The indicator is evaluated at the affine center of each reference triangle.
Faces with ``abs(indicator(center))`` above ``threshold_fraction`` times the
maximum indicator are subdivided. By default, neighboring faces with split
edges are green-refined too, so the adapted reference mesh is conforming. This
mirrors host-mesh adaptation workflows where the high-order curved patches are
rebuilt only after the adapted reference mesh has been created.

Legacy Function
---------------

The older ``integration(...)`` function is still available for notebooks and
existing scripts. New code should prefer ``SurfaceMesh``, ``LevelSetSurface``,
``IntegrationConfig``, and ``integrate`` because these make the numerical setup
explicit and return diagnostic data.
