#distutils: language = c++
#cython: language_level = 3


# # distutils: sources = Rectangle.cpp

import numpy as np
cimport numpy as cnp
from eigency.core cimport *

int_DTYPE   = np.int64
float_DTYPE = np.float

ctypedef cnp.int64_t     int_DTYPE_t
ctypedef cnp.float64_t float_DTYPE_t


from libcpp.vector cimport vector
from libcpp.string cimport string
cnp.import_array()

cdef extern from "FE/Space.h" :
    cdef cppclass Space:
       Space() except +
       int GetNumberOfShapeFunctionsFor(const string&  elemtype);
       void AddDofTo(const string& elemtype, const char& entity, const int& entityNumber, const int& extraKey)
       void Print()
       string  ToStr()

cdef class WrapedSpace:
    cdef Space c_Space
    cdef Space* GetCppObject(self)
