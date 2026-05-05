Projection Diagnostics
======================

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

