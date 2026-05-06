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

``"ModePy_VioreanuRokhlin"``
   Uses the Vioreanu--Rokhlin simplex rule from ModePy. This is the default
   rule. The rule is first mapped from ModePy's biunit triangle to the unit
   reference triangle and is then pulled back through the square-squeezing map.
   In the current ModePy release, this family is available up to degree 20;
   for higher-degree simplex quadrature use ``"ModePy_XiaoGimbutas"``.

``"Pull_back_Gauss"``
   Uses a simplex quadrature rule and pulls the points back through the inverse
   square-squeezing map.

``"Gauss_Legendre"``
   Uses a tensor-product Gauss-Legendre rule on the square.

``"ModePy_XiaoGimbutas"``
   Uses the Xiao--Gimbutas simplex rule from ModePy. The rule is first mapped
   from ModePy's biunit triangle to the unit reference triangle and is then
   pulled back through the square-squeezing map.

``"ModePy_GrundmannMoeller"``
   Uses the Grundmann--Moeller simplex rule from ModePy. For a requested
   integration degree, ``surfgeopy`` chooses the corresponding ModePy order
   whose exactness degree covers the request.

When To Use Which Rule
----------------------

``ModePy_VioreanuRokhlin`` is the package default because it gives a
high-quality simplex cubature family while preserving the square-squeezing
pipeline used by ``surfgeopy``.

``Pull_back_Gauss`` is useful when you want to reproduce the original built-in
simplex rule.

``Gauss_Legendre`` is useful when you want a direct tensor-product rule on the
square parameter domain.

The ``ModePy_*`` rules are useful when you want to compare several established
simplex cubature families while keeping the same cubical reparametrization
pipeline used by the pull-back formulation.

Adding A Rule
-------------

New reference-domain quadrature rules should be added in
``surfgeopy.reference_quadrature`` by returning a ``ReferenceQuadrature`` object.
This keeps the integration kernel independent of the source of the quadrature
points.
