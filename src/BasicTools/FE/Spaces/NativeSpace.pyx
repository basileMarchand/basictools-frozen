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

cimport BasicTools.FE.Spaces.NativeSpace as cNS

cdef class WrapedSpace:
    #cdef Space c_Space hold c++ instance
    cdef Space* GetCppObject(self):
        return &(self.c_Space)

    def SetDataFromPython(self,space):
        for name,data in space.items():
            for nsf in range(data.GetNumberOfShapeFunctions()):
                on,idxI,idxII = data.dofAttachments[nsf]
                if idxI is None :
                    idxI = -1
                if idxII is None :
                    idxII = -1

                #print("name",name)
                #print("on",on)
                #print("idxI",idxI)
                #print("idxII",idxII)
                if on == "F2" :
                    on = "E"
                #print("on",on)
                self.c_Space.AddDofTo(name.encode(),ord(on[0].encode()), int(idxI),int(idxII) )