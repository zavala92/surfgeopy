Quadrature Rules
================

``surfgeopy`` exposes quadrature rules through ``IntegrationConfig``:

.. code-block:: python

   config = IntegrationConfig(
       interpolation_degree=6,
       integration_degree=14,
       quadrature_rule="Gauss_Legendre",
   )

Available Rules
---------------

``"Pull_back_Gauss"``
   Uses a simplex quadrature rule and pulls the points back through the inverse
   square-squeezing map. This is the default rule.

``"Gauss_Legendre"``
   Uses a tensor-product Gauss-Legendre rule on the square.

``"RecursiveNodes_GaussLegendre"``
   Uses the simplex Gauss-Legendre rule from ``recursivenodes`` on the unit
   triangle, then pulls the points back through square-squeezing.

When To Use Which Rule
----------------------

``Pull_back_Gauss`` is a good default for the method implemented in
``surfgeopy``.

``Gauss_Legendre`` is useful when you want a direct tensor-product rule on the
square parameter domain.

``RecursiveNodes_GaussLegendre`` is useful for experimenting with the
``recursivenodes`` quadrature backend while keeping the same high-level
``IntegrationConfig`` interface.

Adding A Rule
-------------

New reference-domain quadrature rules should be added in
``surfgeopy.reference_quadrature`` by returning a ``ReferenceQuadrature`` object.
This keeps the integration kernel independent of the source of the quadrature
points.

