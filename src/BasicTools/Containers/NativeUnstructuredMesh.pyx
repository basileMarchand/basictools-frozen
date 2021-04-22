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

cimport BasicTools.Containers.NativeUnstructuredMesh as cNUM


cdef class CUnstructuredMesh():
    #cdef cNUM.UnstructuredMesh c_UMesh  # Hold a C++ instance which we're wrapping

    cdef cNUM.UnstructuredMesh* GetCppPointer(self) nogil:
        return &(self.cpp_object)

    def SetNodes(self, cnp.ndarray[float_DTYPE_t, ndim=2,mode="c"] nodes not None):
        self.cpp_object.SetNodes(FlattenedMapWithOrder[Matrix, float_DTYPE_t, Dynamic, Dynamic, RowMajor](nodes))

    def SetOriginalIds(self, cnp.ndarray[int_DTYPE_t, ndim=1, mode="c"] originalIds not None):
        self.cpp_object.SetOriginalIds(FlattenedMap[Matrix, int_DTYPE_t, Dynamic, _1](originalIds))

    def AddNodalTag(self, name,  cnp.ndarray[int_DTYPE_t, ndim=1, mode="c"] ids not None):
        self.cpp_object.AddNodalTag(name, FlattenedMap[Matrix, int_DTYPE_t, Dynamic, _1](ids))

    def SetDataFromPython(self,pyUM):
        pyUM.GetPosOfNodes()
        self.SetNodes(pyUM.nodes)
        self.SetOriginalIds(pyUM.originalIDNodes)
        pyUM.nodesTags.Tighten()
        for k in pyUM.nodesTags.keys():
            self.AddNodalTag(k.encode(),pyUM.nodesTags[k].GetIds())

        for k,v in pyUM.elements.items():
            self.cpp_object.AddElemens(k.encode(),
                                    FlattenedMapWithOrder[Matrix, int_DTYPE_t, Dynamic, Dynamic, RowMajor](v.connectivity),
                                    FlattenedMap[Matrix, int_DTYPE_t, Dynamic, _1](v.originalIds) )

            for tn in v.tags.keys():
                self.cpp_object.AddElementTag(k.encode(),tn.encode(),
                                   FlattenedMap[Matrix, int_DTYPE_t, Dynamic, _1](v.tags[tn].GetIds()) )

    def Print(self):
        self.cpp_object.Print()

    def __str__(self):
        res = self.cpp_object.ToStr().decode('UTF-8')
        res += "Done"
        return res

def CheckIntegrity(GUI=False):
    obj = CUnstructuredMesh()
    print(obj)
    return "OK"

if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover
