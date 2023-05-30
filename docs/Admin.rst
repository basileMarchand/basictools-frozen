
.. _fordevs:

**************
For Developers
**************

Testing Infrastructure
######################

BasicTools comes with two ways of executing the automated tests.
Remember that if not all de optional dependencies are installed some tests will fail.
First, pytest [#pytestdoc]_ by simply executing in root directory of the library:

.. code-block::

    pytest

A files ``conftest.py``  and ``pytest.ini`` present at the root of the repository is responsible of the pytest configuration.


Second, a simple in-house tool: every module must have a function called ``CheckIntegrity`` that takes no
argument and returns the string ``"ok"`` if and only if the test was successful.
Any other return value (or a raised exception) will be interpreted as a failed test.

The ``__init__.py`` must have a variable named ``_test`` (the use of the variable ``__all__`` is depreciated) listing all submodules to be tested so that the test infrastructure works as intended.

Some tests need to write data to disk or read data from the test data directory.
Two functions are available to help writing tests :

*  GetTestDataPath() : Returns the path of the data directory (``from BasicTools.TestData import GetTestDataPath`` )
*  TestTempDir: A class to handle the creation of a directory to hold temporary data (``from BasicTools.Helpers.Tests import TestTempDir``)

The function TestAll() (in the module ``BasicTools.Helpers.Tests`` )is used to test the library (see documentation of this function for more information).

This function can be executed using the command:

.. code-block::

    python -m BasicTools.Helpers.Tests

For more in formation about the options use the command :

.. code-block::

    python -m BasicTools.Helpers.Tests -h


Coverage
########

Coverage is an important part of the development process.
To activate the coverage during test use the ``-c`` option.

If you want to tell ``coverage.py`` to ignore some part of the code, use the ``#pragma : no cover`` comment.
See also [#coveragedoc]_.

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
