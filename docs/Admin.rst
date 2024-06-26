
.. _fordevs:

**************
For Developers
**************

Testing Infrastructure
######################

BasicTools comes with two ways of executing the automated tests.
Remember that if not all the optional dependencies are installed some tests will fail.
First, pytest [#pytestdoc]_ by simply executing in root directory of the library:

.. code-block::

    pytest

A files ``conftest.py``  and ``pytest.ini`` present at the root of the repository is responsible of the pytest configuration.

And the in-house testing tool using.

.. code-block::

    python -m BasicTools.Helpers.Tests

This module search all the modules recursively (in BasicTools) and collect all the functions named ``CheckIntegrity``
that takes one boolean argument (True to to generate output in the for of GUI to ease the manual debugging).
The ``CheckIntegrity``  function returns the string ``"ok"`` if and only if the test was successful (or ``"skip"`` to ignore the test).
Any other return value (or a raised exception) will be interpreted as a failed test.

.. note::
    This tool can be used to test a installed version of BasicTools to check the integrity of the installation on the final user environment.


The ``__init__.py`` must have a variable named ``_test`` (the use of the variable ``__all__`` is depreciated) listing all submodules to be tested so that the test infrastructure works as intended.

Some tests need to write data to disk or read data from the test data directory.
Two functions are available to help writing tests :

*  GetTestDataPath() : Returns the path of the data directory (``from BasicTools.TestData import GetTestDataPath`` )
*  TestTempDir: A class to handle the creation of a directory to hold temporary data (``from BasicTools.Helpers.Tests import TestTempDir``)


For more in formation about the options, use the command :

.. code-block::

    python -m BasicTools.Helpers.Tests -h

.. note::
    Some test will fail on some configuration
    Some packages are not available on some platform and some classes must be used inside a specific enviroments (paraview for example)
    The current know issues are :

    -  BasicTools.IO.Catalyst is intended to be used inside a python console in ParaView ([More about](https://www.paraview.org/Wiki/ParaView/Catalyst/Overview)).


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
