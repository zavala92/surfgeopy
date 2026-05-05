Quickstart
==========

This page shows the recommended modern API for integrating a scalar function
over an implicit surface.

The basic workflow is:

1. Load a triangulated reference mesh.
2. Define the implicit surface through ``phi`` and ``grad_phi``.
3. Choose an ``IntegrationConfig``.
4. Call ``integrate`` and read ``result.total``.

Unit Sphere Area
----------------

The unit sphere is represented as the zero level set

.. math::

   \phi(x, y, z) = x^2 + y^2 + z^2 - 1.

Its exact area is :math:`4\pi`.

.. code-block:: python

   import numpy as np

   from surfgeopy import IntegrationConfig, LevelSetSurface, SurfaceMesh, integrate


   def phi(x: np.ndarray) -> float:
       return x[0] ** 2 + x[1] ** 2 + x[2] ** 2 - 1.0


   def grad_phi(x: np.ndarray) -> np.ndarray:
       return np.array([2.0 * x[0], 2.0 * x[1], 2.0 * x[2]])


   mesh = SurfaceMesh.from_mat("tests/mesh_test/sphere_N=104.mat")
   surface = LevelSetSurface(mesh, phi, grad_phi)
   config = IntegrationConfig(
       interpolation_degree=6,
       refinement_level=1,
       integration_degree=14,
       quadrature_rule="Gauss_Legendre",
   )

   result = integrate(surface, lambda _: 1.0, config)
   print(result.total)
   print(result.n_quadrature_points)

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

Legacy Function
---------------

The older ``integration(...)`` function is still available for notebooks and
existing scripts. New code should prefer ``SurfaceMesh``, ``LevelSetSurface``,
``IntegrationConfig``, and ``integrate`` because these make the numerical setup
explicit and return diagnostic data.

