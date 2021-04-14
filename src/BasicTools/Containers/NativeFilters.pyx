#distutils: language = c++
#cython: language_level = 3


import numpy as np
cimport numpy as cnp
from eigency.core cimport *

int_DTYPE   = np.int64
float_DTYPE = np.float

#ctypedef cnp.int64_t     int_DTYPE_t
#ctypedef cnp.float64_t float_DTYPE_t


from libcpp.vector cimport vector
from libcpp.string cimport string
cnp.import_array()


cimport BasicTools.Containers.NativeFilters as cNF

cdef class WrapElementFilterEvaluated:
    #cdef ElementFilterEvaluated cpp_object
    cdef cNF.ElementFilterEvaluated* GetCppObject(self):
        return &(self.cpp_object)

    def callnumpysure(self, elemtype, cnp.ndarray[int_DTYPE_t, ndim=1, mode="c"] ids not None):
        self.cpp_object.SetIdsToTreatFor(elemtype,
                                         FlattenedMap[Matrix, int_DTYPE_t, Dynamic, _1](ids) )

    def SetIdsToTreat(self,mesh,elementFilter ):
        self.GetCppObject().Clear()
        if elementFilter is None:
           for elemtype, data in mesh.elements.items():
               self.callnumpysure(elemtype.encode(),np.arange(data.GetNumberOfElements()) )
        else:
            elementFilter.mesh = mesh
            for elemtype, data in   mesh.elements.items():
                ids= elementFilter.GetIdsToTreat(data)
                ids = np.asarray(ids, dtype=int_DTYPE)
                self.callnumpysure(elemtype.encode(), ids )

    def __str__(self):
        res = self.cpp_object.ToStr().decode('UTF-8')
        res += "Done"
        return res
