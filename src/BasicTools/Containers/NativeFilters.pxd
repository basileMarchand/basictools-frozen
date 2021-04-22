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

cdef extern from "Containers/ElementFilter.h" namespace "BasicTools":
    cdef cppclass ElementFilterBase:
        PlainObjectBase& GetIdsToTreat(const string& elemtype)

    cdef cppclass ElementFilterEvaluated(ElementFilterBase):
       ElementFilterEvaluated() except +
       void SetIdsToTreatFor(string& elemtype, FlattenedMap[Matrix, int_DTYPE_t, Dynamic, _1]& ids )
       void Clear()
       string ToStr() 

cdef class CElementFilterEvaluated:
    cdef ElementFilterEvaluated cpp_object
    cdef ElementFilterEvaluated* GetCppPointer(self) nogil
    