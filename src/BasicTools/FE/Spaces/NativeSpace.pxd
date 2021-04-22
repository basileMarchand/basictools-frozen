#distutils: language = c++
#cython: language_level = 3
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
from libcpp.vector cimport vector
from libcpp.string cimport string

import numpy as np
cimport numpy as cnp
cnp.import_array()

from eigency.core cimport *

from BasicTools.CythonDefs cimport int_DTYPE_t,float_DTYPE_t
from BasicTools.NumpyDefs import int_DTYPE,float_DTYPE

cdef extern from "FE/Space.h" namespace "BasicTools" :
    cdef cppclass Space:
       Space() except +
       int GetNumberOfShapeFunctionsFor(const string&  elemtype);
       void AddDofTo(const string& elemtype, const char& entity, const int& entityNumber, const int& extraKey)
       void Print()
       string  ToStr()

cdef class CSpace:
    cdef Space cpp_object
    cdef Space* GetCppPointer(self)
