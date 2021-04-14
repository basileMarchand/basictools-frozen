#distutils: language = c++
#cython: language_level = 3


# # distutils: sources = Rectangle.cpp

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


cimport BasicTools.Containers.NativeUnstructuredMesh as cNUM


cdef class UnstructuredMesh():
    #cdef cNUM.NativeUnstructuredMesh c_UMesh  # Hold a C++ instance which we're wrapping

    cdef cNUM.NativeUnstructuredMesh* GetCppObject(self):
        return &(self.c_UMesh)

    def SetNodes(self, cnp.ndarray[float_DTYPE_t, ndim=2,mode="c"] nodes not None):
        self.c_UMesh.SetNodes(FlattenedMapWithOrder[Matrix, float_DTYPE_t, Dynamic, Dynamic, RowMajor](nodes))

    def SetOriginalIds(self, cnp.ndarray[int_DTYPE_t, ndim=1, mode="c"] originalIds not None):
        self.c_UMesh.SetOriginalIds(FlattenedMap[Matrix, int_DTYPE_t, Dynamic, _1](originalIds))

    def AddNodalTag(self, name,  cnp.ndarray[int_DTYPE_t, ndim=1, mode="c"] ids not None):
        self.c_UMesh.AddNodalTag(name, FlattenedMap[Matrix, int_DTYPE_t, Dynamic, _1](ids))

    def SetDataFromPython(self,pyUM):
        self.SetNodes(pyUM.nodes)
        self.SetOriginalIds(pyUM.originalIDNodes)
        pyUM.nodesTags.Tighten()
        for k in pyUM.nodesTags.keys():
            self.AddNodalTag(k.encode(),pyUM.nodesTags[k].GetIds())

        for k,v in pyUM.elements.items():
            self.c_UMesh.AddElemens(k.encode(),
                                    FlattenedMapWithOrder[Matrix, int_DTYPE_t, Dynamic, Dynamic, RowMajor](v.connectivity),
                                    FlattenedMap[Matrix, int_DTYPE_t, Dynamic, _1](v.originalIds) )

            for tn in v.tags.keys():
                self.c_UMesh.AddElementTag(k.encode(),tn.encode(),
                                   FlattenedMap[Matrix, int_DTYPE_t, Dynamic, _1](v.tags[tn].GetIds()) )

    def Print(self):
        self.c_UMesh.Print()

    def __str__(self):
        res = self.c_UMesh.ToStr().decode('UTF-8')
        res += "Done"
        return res

def CheckIntegrity(GUI=False):
    obj = UnstructuredMesh()
    print(obj)
    return "OK"

if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover
