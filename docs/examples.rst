Examples
========

The examples are grouped by their role in the numerical story. They focus on
the stable smooth surface-quadrature workflow.

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

``examples/adapted_mesh_Gauss_Bonnet_benchmark.ipynb``
   Indicator-refined Gauss-Bonnet benchmark on a biconcave surface.
