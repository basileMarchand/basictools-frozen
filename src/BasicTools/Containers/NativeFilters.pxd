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
#cdef extern from "Containers/ElementFilter.cpp" :
#    pass

cdef extern from "Containers/ElementFilter.h" :
    cdef cppclass ElementFilterBase:
        PlainObjectBase& GetIdsToTreat(const string& elemtype)

       #ElementFilterEvaluated() except +
       #void SetIdsToTreatFor(string& elemtype, FlattenedMap[Matrix, int_DTYPE_t, Dynamic, _1]& ids)

    cdef cppclass ElementFilterEvaluated(ElementFilterBase):
       ElementFilterEvaluated() except +
       void SetIdsToTreatFor(string& elemtype, FlattenedMap[Matrix, int_DTYPE_t, Dynamic, _1]& ids )

       void Clear()
       string ToStr()

cdef class WrapElementFilterEvaluated:
    cdef ElementFilterEvaluated cpp_object
    cdef ElementFilterEvaluated* GetCppObject(self)
    #def SetIdsToTreat(self,mesh,elementFilter )