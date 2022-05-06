***************
Administration
***************

Repositories
############

An open source version of BasicTools is released at regular intervals (around every 3 months) on Gitlab.com [#gitlaburlpublic]_.

Installation
############

Requirements
************

Python minimal version: 3.7


Python packages:

* numpy
* scipy
* sympy
* pyparsing
* cython
* scikit-learn
* scikit-sparse
* vtk
* eigency=1.78
* mkl
* mkl-include
* psutil

C++ packages:

* Eigen [#eigenurl]_

Optionals Python packages (some functionalities may not be available):

* pyamg
* h5py
* meshio
* sphinx
* sphinx-rtd-theme
* setuptools-scm
* pyvista


Optionals Proprietary packages (some functionalities may not be available):

* odbAccess
* abaqusConstants

Conda Installation
******************

The Anaconda [#anacondaurl]_  python package manager provides a, ready to use, package of BasicTools.
The BasicTools Anaconda package [#basictoolsanaconda]_ is managed by the community.

The installation can be done using the following command.


    ``> conda install -c conda-forge basictools``

Manual Installation
*******************

for manual installation:

    > pip install .

In the case the user want to make changes to BasicTools, an installation from sources is mandatory.
The sources can be downloaded from either one of the gitlab repositories.
Then inside the repository folder the user must compile the c++ extensions to take profit of optimized algorithms.

    > python setup.py build_clib
    > python setup.py build_ext --inplace

Then the user is responsible to add the ``BASICTOOLS_REPOSITORY/src/`` folder to the ``PYTHONPATH`` environment variables (more information on [#pythonpathdoc]_).

Or using pip for development:

    > pip install -e .


For Windows
^^^^^^^^^^^

Scikit-sparse package not available in anaconda for windows, should be able to compile it for windows following [#scikitwindows]_.
I good stating point for the installation is https://github.com/EmJay276/scikit-sparse .

Extra pre-requirement:
    - Microsoft Build Tools for C++.
    - suitesparse package

To contribute to BasicTools
###########################

To contribute to BasicTools you can fill a bug report on either of the gitlab with a *minimal non working example*.
This makes the debugging in ours side easier.

If you want to contribute with code you must.

*  clone the master branch of BasicTools from either one of the two gitlab repositories
*  make a development branch
*  modify/created changes, commit changes
*  compile BasicTools
*  test you branch (see section `Testing Infrastructure`_ )
*  accept the Contribution Agreement (see section `Licensing and External Contributions`_ )
*  push your branch to the server
*  create a merge request (on the web)

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

Licensing and External Contributions
####################################

BasicTools license.

.. literalinclude:: ../LICENSE.txt

For contribution the user must accept the Contribution Agreement. This means the user transfer the property of the contributed code to Safran.

.. literalinclude:: ../CONTRIBUTING.md

.. rubric:: Footnotes
.. [#gitlaburlpublic]  https://gitlab.com/drti/basic-tools
.. [#anacondaurl] https://anaconda.org/
.. [#basictoolsanaconda] https://anaconda.org/conda-forge/basictools
.. [#pythonpathdoc] https://docs.python.org/3/using/cmdline.html\\#envvar-PYTHONPATH
.. [#eigenurl] http://eigen.tuxfamily.org
.. [#pytestdoc] https://docs.pytest.org/en/stable/
.. [#coveragedoc] https://coverage.readthedocs.io/en/coverage-5.5/excluding.html
.. [#scikitwindows] https://github.com/xmlyqing00/Cholmod-Scikit-Sparse-Windows