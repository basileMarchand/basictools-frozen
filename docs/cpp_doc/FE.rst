********************
Finite Element Tools
********************

DofNumbering
------------

.. doxygenclass:: BasicTools::DofNumbering

Integration Rules
-----------------

.. doxygenclass:: BasicTools::IntegrationRule
.. doxygenclass:: BasicTools::SpaceIntegrationRule

.. doxygenvariable:: BasicTools::IntegrationRulesAlmanac
.. doxygenfunction:: BasicTools::GetPythonDefinedIntegrationRules

Finite Element Spaces
---------------------

.. doxygenfunction:: BasicTools::GetSpaceFor
.. doxygenfunction:: BasicTools::GetAvailableSpaces

Integration
-----------

.. doxygenfunction:: BasicTools::solve(T&, MapMatrixDDD&)
.. doxygenfunction:: BasicTools::solve(MatrixDDD&, MapMatrixDDD&)


.. doxygenstruct:: BasicTools::LocalSpace
.. doxygenstruct:: BasicTools::MonoElementsIntegralCpp

Weak form classes
-----------------

.. doxygenstruct:: BasicTools::WeakTerm
.. doxygenstruct:: BasicTools::WeakMonom
.. doxygenstruct:: BasicTools::WeakForm

Spaces
------

.. doxygenstruct:: BasicTools::ElementSpace
.. doxygenclass:: BasicTools::Space

.. doxygenclass:: BasicTools::SpaceAtIP

.. doxygenfunction:: BasicTools::EvaluateSpaceAt(const Space &, const SpaceIntegrationRule& )
.. doxygenfunction:: BasicTools::EvaluateSpaceAt(const ElementSpace&, const IntegrationRule& )