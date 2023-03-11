
.. _fordevs:

**************
For Developers
**************

Testing Infrastructure
######################

BasicTools came with a very basic but useful testing infrastructure.
To execute all test you can use either pytest [#pytestdoc]_ or the in house infra.
A file ``conftest.py`` present at the root of the repository is responsible of the pytest configuration.

Every module must have a function called ``CheckIntegrity`` that takes no argument and returns the string ``"ok"`` if and only if the test was successful.
Any other return value (or a raised exception) will be interpreted as a failed test.

The ``__init__.py`` must have a variable named ``_test`` (the use of the variable ``__all__`` is depreciated) listing all submodules to be tested so that the test infrastructure works as intended.

Some tests need to write data to disk or read data from the test data directory.
Two functions are available to help writing tests :

*  GetTestDataPath() : Returns the path of the data directory (``from BasicTools.TestData import GetTestDataPath`` )
*  TestTempDir: A class to handle the creation of a directory to hold temporary data (``from BasicTools.Helpers.Test import TestTempDir``)

The function TestAll() (in the module ``BasicTools.Helpers.Test`` )is used to test the library (see documentation of this function for more information).

This function can be executed using the command:

    ``> python -m BasicTools.Helpers.Test``

For more in formation about the options use the command :

    ``> python -m BasicTools.Helpers.Test -h``


Coverage
########

Coverage is an important part of the development process.
To activate the coverage during test use the ``-c`` option.

If you want to tell ``coverage.py`` to ignore some part of the code, use the ``#pragma : no cover`` comment.
See also : [#coveragedoc]_.

Disabling Tests
###############

Some tests can be disabled using an environment variable.
A typical use case arises when a test relies on an external dependency that may not be available.

The feature relies on the definition of non-empty enviromennt variables in the form : ``"appsname_NO_FAIL"``.

An example is available in the file ``BasicTools/FE/ZmatFemProblem.py`` at the beginning of the ``CheckIntegrity()`` function.

.. rubric:: Footnotes
.. [#basictoolsanaconda] https://anaconda.org/conda-forge/basictools
.. [#pytestdoc] https://docs.pytest.org/en/stable/
.. [#coveragedoc] https://coverage.readthedocs.io/en/coverage-5.5/excluding.html
.. [#scikitwindows] https://github.com/xmlyqing00/Cholmod-Scikit-Sparse-Windows