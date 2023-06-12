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

from BasicTools.CythonDefs cimport CBasicIndexType, CBasicFloatType


cimport BasicTools.Containers.NativeUnstructuredMesh as cNUM
cimport BasicTools.Containers.NativeFilters as cNF
cimport BasicTools.FE.Spaces.NativeSpace as cNS

import BasicTools.Containers.NativeUnstructuredMesh as NUM
import BasicTools.FE.Spaces.NativeSpace as NS

cdef extern from "FE/DofNumbering.h" namespace "BasicTools" :
    cdef cppclass DofNumbering:
        DofNumbering() except +
        CBasicIndexType GetSize()
        CBasicIndexType GetDofOfPoint(int)
        int GetFromConnectivity()
        void SetFromConnectivity(bool )
        void ComputeNumberingFromConnectivity(cNUM.UnstructuredMesh& ) nogil
        void ComputeNumberingGeneral(cNUM.UnstructuredMesh&, cNS.Space&, cNF.ElementFilterBase&) nogil
        void computeDofToPoint()

        MatrixXi & GetNumberingFor(const string & elemtype)
        CBasicIndexType GetSizeOfDofToPoint()
        PlainObjectBase & GetdoftopointLeft()
        PlainObjectBase & GetdoftopointRight()

        void computeDofToCell(cNUM.UnstructuredMesh)
        CBasicIndexType GetSizeOfDofToCell()
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

    cdef DofNumbering* GetCppPointer(self) nogil


