
1) NOTES FOR ALL PEOPLE

Please configure your git before doing anything!!!!!!

  git config --global user.name "Billy Everyteen"
  git config --global user.email "your_email@example.com"

To clone this repository :

  git -c http.sslVerify=false clone https://dXXXXXX@sc-crt-1.safran/git-ms/BasicTools
  XXXXX to be replaced with your safran number

To configure the tools :

For Python : Add the python directory to your PYTHONPATH

The page for the bug repports and more:

	https://sc-crt-1.safran/redmine/
	Project OTTools

2) DEPENDENCIES

    Numpy
    Scipy
    Qt
    Vtk
    Cython


3) INSTALLATION

    to compile the c++ sources run command in the root directory

     > python setup.py build_ext --inplace


4) NOTES FOR CONTRIBUTORS

    We recomend to use the "simple"  behavior for pushing :

       git config --global push.default simple


    Please Read the following page before contributing code :

        https://www.python.org/dev/peps/pep-0008/

    For the moment the only differences from PEP 0008 are :

      Function Names: "CamelCase"

        CamelCase starting with uppercase

      Variables Names: "camelCase"

        camelCase starting with lowercase

    For the moment you must bypass the ssl certificate verification, use the
    following command to push your changes:

        git -c http.sslVerify=false push


5) TESTING INFRASTRUCTURE

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




