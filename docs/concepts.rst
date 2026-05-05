Core Concepts
=============

``surfgeopy`` computes high-order surface integrals on smooth embedded
surfaces when an implicit representation is known.

Reference Mesh
--------------

The input mesh is a triangulated reference geometry. It does not need to be a
high-order mesh. Each reference triangle is used as the starting point for a
curved approximation of the true surface.

Implicit Surface
----------------

The surface is represented by a level-set function and its gradient:

.. math::

   S = \{x \in \mathbb{R}^3 : \phi(x) = 0\}.

The gradient is used by the closest-point projection routine to move
interpolation nodes from the reference mesh onto the implicit surface.

Curved Geometry Approximation
-----------------------------

For each reference triangle, ``surfgeopy`` maps interpolation nodes through the
square-squeezing transform and projects them to the implicit surface. These
projected nodes define a polynomial approximation of the curved surface patch.

Quadrature
----------

After the curved patch is built, the integrand and geometric Jacobian are
evaluated at quadrature points. The final integral is the weighted sum over all
patches.

What Information Is Required?
-----------------------------

High-order accuracy requires curved-surface information. In ``surfgeopy`` that
information comes from the implicit representation: ``phi`` and ``grad_phi``.
A triangle mesh alone only describes a piecewise-linear surface and cannot, by
itself, recover the same geometry accuracy.

