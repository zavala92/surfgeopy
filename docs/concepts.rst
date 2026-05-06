Core Concepts
=============

``surfgeopy`` computes high-order surface integrals on smooth embedded
surfaces when an implicit representation is available. The package starts from
a linear triangulated approximation and constructs curved surface patches by
interpolating the closest-point projection. The interpolation task is pulled
back from each triangle to the reference square, where tensor-product
Chebyshev--Lobatto grids provide stable high-order interpolation nodes.

Reference Mesh
--------------

The input mesh is a triangulated reference geometry. It does not need to be a
high-order mesh: each triangle is used as a geometric starting point for
constructing a curved patch on the target surface.

For a triangle :math:`T_i` in the reference triangulation :math:`T`, let

.. math::

   \tau_i : \Delta_2 \rightarrow T_i

be the affine map from the reference triangle :math:`\Delta_2` to the physical
mesh triangle.

Implicit Surface
----------------

The target surface is represented by a level-set function and its gradient:

.. math::

   S = \{x \in \mathbb{R}^3 : \phi(x) = 0\}.

The gradient is used by the closest-point projection routine. Given a point on
or near the linear reference mesh, the projection moves it onto the implicit
surface. For the patch associated with :math:`T_i`, denote this projection by

.. math::

   \pi_i : T_i \rightarrow S_i.

Square-Squeezing Pullback
-------------------------

The defining step of the method is the pullback of the surface interpolation
problem to the square. Let

.. math::

   \sigma : \square_2 \rightarrow \Delta_2

be the square-squeezing map from the reference square to the reference
triangle. ``surfgeopy`` represents the curved surface patch by the composed map

.. math::

   \varphi_i : \square_2 \rightarrow S_i,
   \qquad
   \varphi_i = \pi_i \circ \tau_i \circ \sigma.

Thus, instead of interpolating directly on a triangle, the package interpolates
the geometry map :math:`\varphi_i` on the tensor-product domain
:math:`\square_2`.

Chebyshev--Lobatto Geometry Interpolation
-----------------------------------------

On the square, ``surfgeopy`` uses Chebyshev--Lobatto grids to construct a
vector-valued tensor-polynomial interpolant of the surface map:

.. math::

   Q_{G_{2,k}} \varphi_i
   = \sum_{\alpha \in A_{2,k}} b_\alpha N_\alpha.

Here :math:`Q_{G_{2,k}} \varphi_i` is the degree-:math:`k` interpolant on the
grid :math:`G_{2,k}`, and the coefficients
:math:`b_\alpha \in \mathbb{R}` of the Newton interpolation can be computed in
closed form. Chebyshev--Lobatto nodes are used because they provide stable
high-order interpolation and help avoid the Runge phenomenon.

Surface Integral Approximation
------------------------------

Replacing :math:`\varphi_i` by :math:`Q_{G_{2,k}}\varphi_i` yields the
high-order geometric approximation used in the surface integral. For an
integrand :math:`f`, ``surfgeopy`` evaluates

.. math::

   \int_S f\,dS
   \approx
   \sum_{i=1}^K
   \int_{\square_2}
   (f \circ \varphi_i)(\mathrm{x})
   \sqrt{
      \det\left(
         (DQ_{G_{2,k}}\varphi_i(\mathrm{x}))^T
         DQ_{G_{2,k}}\varphi_i(\mathrm{x})
      \right)
   }
   d\mathrm{x}.

The remaining integral is computed by high-order quadrature:

.. math::

   \sum_{i=1}^K
   \sum_{\mathrm{p} \in P}
   \omega_{\mathrm{p}}
   (f \circ \varphi_i)(\mathrm{p})
   \sqrt{
      \det\left(
         (DQ_{G_{2,k}}\varphi_i(\mathrm{p}))^T
         DQ_{G_{2,k}}\varphi_i(\mathrm{p})
      \right)
   }.

Quadrature Choices
------------------

After the curved patch has been built, ``surfgeopy`` evaluates the integrand
and geometric Jacobian at quadrature points. The package supports two natural
families of quadrature rules:

``Gauss_Legendre``
   A tensor-product Gauss--Legendre rule on the square
   :math:`\square_2`.

``Pull_back_Gauss``
   A simplex quadrature rule on :math:`\Delta_2`, pulled back to
   :math:`\square_2` using the inverse of the square-squeezing map. This keeps
   the rule aligned with the original triangulation while still using the
   square-based interpolation representation.

Why the Implicit Representation Matters
---------------------------------------

High-order accuracy requires high-order geometric information. In
``surfgeopy`` that information is supplied by the implicit representation:
``phi`` and ``grad_phi`` define the closest-point projection used to recover
curved surface patches from the linear mesh. A triangle mesh alone describes
only a piecewise-linear surface and cannot, by itself, recover the same
high-order geometry accuracy.
