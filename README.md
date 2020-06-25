1) DEPENDENCIES

    python minimal version: 3.6

    PYTHON OPEN-SOURCE DEPENDENCIES
    numpy
    scipy
    scikit-sparse
    matplotlib
    vtk
    sympy
    pyamg
    h5py
    pyparsing
    Cython
    sphinx
    pytest
    pytest-cov
    setuptools

    C++ OPEN-SOURCE DEPENDENCIES
    Eigen (http://eigen.tuxfamily.org)

    THIRD-PARTY PROPRIETARY DEPENDENCIES
    odbAccess
    abaqusConstants

    FOR WINDOWS:
      install Microsoft Visual C++ Build Tools to use eigen,
      scikit-sparse package not available in anaconda for windows, should be able to compile it for windows following https://github.com/xmlyqing00/Cholmod-Scikit-Sparse-Windows,
      tested on windows10: 99 tests OK, 3 tests not OK (Linalg.LinearSolver, IO.Wormhole (SIGALARM not supported on windows), IO.CodeInterface



2) INSTALLATION

    SETUP:

    An environement variable with the path to the EIGEN library must be defined:

     > export EIGEN_INC=/Path/To/Eigen/Library

    COMPILATION:,

    Run the following command in the root directory :

     > python setup.py build_ext --inplace


3) TESTING INFRASTRUCTURE

    Every module must have a function called "CheckIntegrity" that takes no
    argument and returns the string "ok" if and only if the test was successful.

    The __init__.py must have a variable named _test (the use of the variable
    __all__ is depreciated) listing all submodules to be tested so that the test
    infrastructure works as intended.

    Two functions are available to help writing tests :

    -   GetTestDataPath() : Returns the path of the data directory
    -   TestTempDir(): Returns a directory to hold temporary data

    The function TestAll() is used to test the library (see documentation of
    this function for more information).

    COVERAGE :

    If you want to tell coverage.py to ignore some part of the code, use the
    "#pragma : no cover" comment. See also :
    http://coverage.readthedocs.org/en/coverage-4.0.3/excluding.html

    DISABLING TESTS :

    Some tests can be disabled using an environment variable. A typical use
    case arises when a test relies on an external dependency that may not be
    available.

    The feature relies on the definition of non-empty enviromennt variables of
    the form :

        "appsname_NO_FAIL"

    An example is available in the file BasicTools/FE/ZmatFemProblem.py at the
    beginning of the CheckIntegrity() function.
