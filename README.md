1) DEPENDENCIES

    python minimal version: 3.7

    PYTHON OPEN-SOURCE DEPENDENCIES
    numpy
    scipy
    sympy
    pyparsing
    cython
    scikit-learn
    scikit-sparse
    vtk
    eigency=1.78
    mkl
    mkl-include
    psutil
    pyamg
    h5py
    meshio
    sphinx
    sphinx-rtd-theme
    setuptools-scm
    pyvista
    psutil
    networkx
    pywin32 [Only for windows]


    C++ OPEN-SOURCE DEPENDENCIES
    Eigen (http://eigen.tuxfamily.org)

    THIRD-PARTY PROPRIETARY DEPENDENCIES
    odbAccess
    abaqusConstants

    FOR WINDOWS:
      install Microsoft Visual C++ Build Tools to use eigen,
      scikit-sparse package not available in anaconda for windows, some functionality will be missing,

2) INSTALLATION

    FOR USER:

    For conda you can create packages using the recipes available on the sources:

        recipes/compiled/  -> for compiled version of BasicTools
        recipes/noarch/  -> (not recommended) for a pure python (slower) version of BasicTools

    FOR DEVELOPERS:

    For development using a conda or other type of environment manager :

        create an environment with all the requirements


    if the compilation script cant find eigen please add the environment  variable with the path to the EIGEN library must be defined (normally you don't have to do this):

        > export EIGEN_INC=/Path/To/Eigen/Library

    COMPILATION:

    Run the following command in the root directory :

        > python setup.py build_clib
        > python -m pip install --no-deps  -e . -vv

    This will install the library in developer mode (-e), to reinstall please remove the "build" before running
    the commands again.

    WARNING: if you use a shared conda environnement, then you must only compile (not install) BasicTools

        > python setup.py build_clib
        > python setup.py build_ext --inplace
        > export PYTHONPATH=${PYTHONPATH}:/path/to/BasicTools/src/
        (for bash)


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

    With a correct configuration of the environment, the following command should
    return 100% of successful tests:

        > python -m BasicTools -k Catalyst

    COVERAGE :

    If you want to tell coverage.py to ignore some part of the code, use the
    "#pragma : no cover" comment. See also :
    http://coverage.readthedocs.org/en/coverage-4.0.3/excluding.html

    DISABLING TESTS :

    Some tests can be disabled using an environment variable. A typical use
    case arises when a test relies on an external dependency that may not be
    available.

    The feature relies on the definition of non-empty environment variables of
    the form :

        "appsname_NO_FAIL"

    An example is available in the file BasicTools/FE/ZmatFemProblem.py at the
    beginning of the CheckIntegrity() function.

4) DOCUMENTATION

    The documentation for BasicTools can be compiled using sphinx

        > python setup.py build_sphinx

    Also the documentation can be found at

        https://basictools.readthedocs.io/en/latest/