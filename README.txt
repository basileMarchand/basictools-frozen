1) DEPENDENCIES

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
    coverage

    C++ OPEN-SOURCE DEPENDENCIES
    Eigen (http://eigen.tuxfamily.org)

    THIRD-PARTY PROPRIETARY DEPENDENCIES
    odbAccess
    abaqusConstants

2) INSTALLATION

    SETUP:

    An environement variable with the path to the EIGEN library must be defined:

     > export EIGEN_INC=/Path/To/Eigen/Library

    COMPILATION:,

    Run the following command in the root directory :

     > python setup.py build_ext --inplace

3) NOTES FOR CONTRIBUTORS

    Please read the following page before contributing code :

        https://www.python.org/dev/peps/pep-0008/

    The only deviations from PEP 8 are the following :

      Function Names: "CamelCase"

        CamelCase starting with uppercase

      Variables Names: "camelCase"

        camelCase starting with lowercase

4) TESTING INFRASTRUCTURE

    Every module must have a function called "CheckIntegrity" that takes no
    argument and returns the string "ok" if and only if the test was successful.

    The __init__.py must have a variable named __all__ listing all submodules
    so that the test infrastructure works as intended.

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
