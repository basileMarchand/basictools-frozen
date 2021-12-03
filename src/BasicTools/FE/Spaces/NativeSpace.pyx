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

from BasicTools.CythonDefs cimport CBasicIndexType, CBasicFloatType
from BasicTools.NumpyDefs import PBasicIndexType, PBasicFloatType

cimport BasicTools.FE.Spaces.NativeSpace as cNS

cdef class CSpace:

    cdef cNS.Space* GetCppPointer(self):
        return &(self.cpp_object)

    def SetDataFromPython(self,space):
        for name,data in space.items():
            for nsf in range(data.GetNumberOfShapeFunctions()):
                on,idxI,idxII = data.dofAttachments[nsf]
                if idxI is None :
                    idxI = -1
                if idxII is None :
                    idxII = -1
                if on == "F2" :
                    on = "E"
                self.cpp_object.AddDofTo(name.encode(),ord(on[0].encode()), PBasicIndexType(idxI), PBasicIndexType(idxII) )
