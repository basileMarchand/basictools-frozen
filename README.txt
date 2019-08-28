1) DEPENDENCIES

    Numpy
    Scipy
    Qt
    Vtk
    Cython
    Eigen

2) INSTALLATION

    to compile the c++ sources run command in the parent directory of BasicTools

     > python BasicTools/setup.py build_ext --inplace

3) NOTES FOR CONTRIBUTORS

    Please Read the following page before contributing code :

        https://www.python.org/dev/peps/pep-0008/

    For the moment the only differences from PEP 0008 are :

      Function Names: "CamelCase"

        CamelCase starting with uppercase

      Variables Names: "camelCase"

        camelCase starting with lowercase

4) TESTING INFRASTRUCTURE

    Every module must have a function called "CheckIntegrity", this function has
    no arguments and must return the string "OK" if the test was successful.

    The __init__.py must have a variable named __all__ contanig the list of all
    the submodules for the test infrastructure to work.

    Two functions are avilable to help writing test:

    -   GetTestDataPath() : Function to get the path of the data directory
    -   TestTempDir(): Function to get a temporary directory (to store temp data)

    To test the library the function TestAll() is used (see doc of this function
    for more information).

    COVERAGE:  If you want to tell coverage.py to ignore some part of the code,
               use the "#pragma : no cover" comment. More information about
               coverage http://coverage.readthedocs.org/en/coverage-4.0.3/excluding.html

    DISABLING TESTS :

    Some tests can be disabled using an enviroment variable. For Example in the case
    the test needs an external program not available in the current system.

    The idea is to use a non empty enviroment variable in the form of :

        "appsname_NO_FAIL"

    Please look at the file BasicTools/FE/ZmatFemProblem.py at the beginning of
    the CheckIntegrity() function for an example.
