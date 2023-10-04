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

from BasicTools.CythonDefs cimport  CBasicIndexType, CBasicFloatType
from BasicTools.NumpyDefs import int_DTYPE,float_DTYPE

cdef extern from "Containers/UnstructuredMesh.h" namespace "BasicTools":
    cdef cppclass UnstructuredMesh:
       UnstructuredMesh() except +
       void SetOriginalIds(FlattenedMap[Matrix, CBasicIndexType, Dynamic, _1]& )
       void SetNodes(FlattenedMapWithOrder[Matrix, CBasicFloatType, Dynamic, Dynamic, RowMajor]& )
       void AddNodalTag(string& name, FlattenedMap[Matrix, CBasicIndexType, Dynamic, _1] &arg1)
       void AddElements(string& name,
                       FlattenedMapWithOrder[Matrix, CBasicIndexType, Dynamic, Dynamic, RowMajor] &conn ,
                       FlattenedMap[Matrix, CBasicIndexType, Dynamic, _1] &ids)

       void AddElementTag(string& elementType,
                          string& tagname,
                          FlattenedMap[Matrix, CBasicIndexType, Dynamic, _1] &ids)

       void Print()
       string  ToStr()


cdef class CUnstructuredMesh:
    cdef UnstructuredMesh cpp_object
    cdef UnstructuredMesh* GetCppPointer(self) nogil
