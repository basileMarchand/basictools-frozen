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

cdef extern from "Containers/UnstructuredMesh.h" :
    cdef cppclass NativeUnstructuredMesh:
       NativeUnstructuredMesh() except +
       void SetOriginalIds(FlattenedMap[Matrix, int_DTYPE_t, Dynamic, _1]& )
       void SetNodes(FlattenedMapWithOrder[Matrix, float_DTYPE_t, Dynamic, Dynamic, RowMajor]& )
       void AddNodalTag(string& name, FlattenedMap[Matrix, int_DTYPE_t, Dynamic, _1] &arg1)
       void AddElemens(string& name,
                       FlattenedMapWithOrder[Matrix, int_DTYPE_t, Dynamic, Dynamic, RowMajor] &conn ,
                       FlattenedMap[Matrix, int_DTYPE_t, Dynamic, _1] &ids)

       void AddElementTag(string& elementType,
                          string& tagname,
                          FlattenedMap[Matrix, int_DTYPE_t, Dynamic, _1] &ids)

       void Print()
       string  ToStr()


cdef class UnstructuredMesh:
    cdef NativeUnstructuredMesh c_UMesh
    cdef NativeUnstructuredMesh* GetCppObject(self)
