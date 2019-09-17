Please read the following page before contributing code :

    https://www.python.org/dev/peps/pep-0008/

The only deviations from PEP 8 are the following :

    Function Names: "CamelCase"

    CamelCase starting with uppercase

    Variables Names: "camelCase"

    camelCase starting with lowercase

    Line limit length of 79 characters


The following rules have been followed:

a). API are kept simple and evolutive.

b). Data structures are defined in Containers, we limit the number of advanced 
functions in these classes at maximum. This advanced functions are defined
in other files.

c). Unitary tests are in Containers classes and functional tests in files
containing more advanced classes. The latter also serve as exemples and
client code for the library.

d). Each file should at most one class (even if there are already some
deviations...).

e). CheckIntegrities must be local (as litte imports as possible), aiming to test
only the functions defined in the file, as much a spossible. All functions
in the file must be tested, if possible. Please limit the use of functions
from other file to reach that goal, and use small and simple data, otherwise
changes will be painful to propagate if many CheckIntegrities must be updated
as well.

f). Coverage must be kept higher than 80%. 

g). Favor imports at the beginning of files (to be available for CheckIntegrities
as well without repeat).
    
