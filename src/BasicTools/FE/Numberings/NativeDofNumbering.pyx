#distutils: language = c++
#cython: language_level = 3
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool

import numpy as np
cimport numpy as cnp
cnp.import_array()

from eigency.core cimport *

from BasicTools.CythonDefs cimport int_DTYPE_t,float_DTYPE_t
from BasicTools.NumpyDefs import int_DTYPE,float_DTYPE


cimport BasicTools.Containers.NativeUnstructuredMesh as cNUM
cimport BasicTools.Containers.NativeFilters as cNF
cimport BasicTools.FE.Spaces.NativeSpace as cNS

import BasicTools.Containers.NativeUnstructuredMesh as NUM
import BasicTools.FE.Spaces.NativeSpace as NS

cdef extern from "FE/DofNumbering.h" namespace "BasicTools" :
    cdef cppclass DofNumbering:
        DofNumbering() except +
        int GetSize()
        int GetDofOfPoint(int)
        int GetFromConnectivity()
        void SetFromConnectivity(bool )
        void ComputeNumberingFromConnectivity(cNUM.UnstructuredMesh& )
        void ComputeNumberingGeneral(cNUM.UnstructuredMesh&, cNS.Space&, cNF.ElementFilterBase&)
        void computeDofToPoint()

        PlainObjectBase & GetNumberingFor(const string & elemtype) 
        int_DTYPE_t GetSizeOfDofToPoint()
        PlainObjectBase & GetdoftopointLeft()
        PlainObjectBase & GetdoftopointRight()

        void computeDofToCell(cNUM.UnstructuredMesh)
        int_DTYPE_t GetSizeOfDofToCell()
        PlainObjectBase & GetdoftocellLeft()
        PlainObjectBase & GetdoftocellRight()
        bool HasNumberingFor(const string & elemtype)
        string ToStr()

cdef class NativeDofNumbering:
    cdef DofNumbering dn_cpp
    cdef cnp.ndarray _doftopointLeft
    cdef cnp.ndarray _doftopointRight
    cdef cnp.ndarray _doftocellLeft
    cdef cnp.ndarray _doftocellRight
    cdef bool doftopointComputed
    cdef dict __dict__
    #def get(self,key,default=None):
    #    return self.numbering.get(key,default)
    def __init__(self):
        self.doftopointComputed = False
        self.mesh = None
    def __getitem__(self,key):
        if key == "size":
           # print("Please use the new API of DofNumbering : DofNumbering.size")
            return self.size

        if key == "fromConnectivity":
           # print("Please use the new API of DofNumbering : DofNumbering.fromConnectivity")
            return self.fromConnectivity

        if key == "doftopointLeft":
            return self.doftopointLeft

        if key == "doftopointRight":
            return self.doftopointRight

        if key == "doftocellLeft":
            return self.doftocellLeft

        if key == "doftocellRight":
            return self.doftocellRight


        if self.dn_cpp.HasNumberingFor(key.encode()):
            return  ndarray_view(self.dn_cpp.GetNumberingFor(key.encode()))
        else:
            return None
            #raise KeyError

    def get(self, key, default):
        return self[key]

    @property
    def size(self):
        return self.dn_cpp.GetSize()

    @property
    def fromConnectivity(self):
        return self.dn_cpp.GetFromConnectivity()

    def GetDofOfPoint(self,val):
        return self.dn_cpp.GetDofOfPoint(val)

    #@cython.boundscheck(False)  # Deactivate bounds checking
    #@cython.wraparound(False)   # Deactivate negative indexing.
    def ComputeNumberingFromConnectivity(self,mesh,notused):
        cdef cNUM.CUnstructuredMesh obj = NUM.CUnstructuredMesh()
        obj.SetDataFromPython(mesh)
        self.dn_cpp.ComputeNumberingFromConnectivity(obj.GetCppPointer()[0])
        return self

    def ComputeNumberingGeneral(self,mesh,space,elementFilter=None,discontinuous=False):
        #self.fromConnectivity = False
        self.mesh = mesh

        cdef cNUM.CUnstructuredMesh obj = NUM.CUnstructuredMesh()
        obj.SetDataFromPython(mesh)

        cdef cNS.CSpace s = cNS.CSpace()
        s.SetDataFromPython(space)

        cdef cNF.CElementFilterEvaluated ef = cNF.CElementFilterEvaluated()
        ef.SetIdsToTreat(mesh,elementFilter)
        self.dn_cpp.ComputeNumberingGeneral(obj.GetCppPointer()[0], s.GetCppPointer()[0], ef.GetCppPointer()[0])

        return self

    @property
    def doftopointLeft(self):
        if self.dn_cpp.GetSizeOfDofToPoint() == 0:
            return np.array([],dtype=int)
        else:
            return  ndarray_view(self.dn_cpp.GetdoftopointLeft())

    @property
    def doftopointRight(self):
        if self.dn_cpp.GetSizeOfDofToPoint() == 0:
            return np.array([],dtype=int)
        else:
            return  ndarray_view(self.dn_cpp.GetdoftopointRight())

    @property
    def doftocellLeft(self):
        cdef cNUM.CUnstructuredMesh obj = NUM.CUnstructuredMesh()
        obj.SetDataFromPython(self.mesh)
        self.dn_cpp.computeDofToCell(obj.GetCppPointer()[0] )
        if self.dn_cpp.GetSizeOfDofToCell() == 0:
            return np.array([],dtype=int)
        else:
            return ndarray(self.dn_cpp.GetdoftocellLeft())

    @property
    def doftocellRight(self):
        cdef cNUM.CUnstructuredMesh obj = NUM.CUnstructuredMesh()
        obj.SetDataFromPython(self.mesh)
        self.dn_cpp.computeDofToCell(obj.GetCppPointer()[0] )
        if self.dn_cpp.GetSizeOfDofToCell() == 0:
            return np.array([],dtype=int)
        else:
            return  ndarray_view(self.dn_cpp.GetdoftocellRight())

    def __str__(self):
        res = self.dn_cpp.ToStr().decode('UTF-8')
        return res

def CheckIntegrity(GUI=False):
    import BasicTools.FE.DofNumbering  as DN
    return DN.CheckIntegrityUsingAlgo("DictBase",GUI)

if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover
