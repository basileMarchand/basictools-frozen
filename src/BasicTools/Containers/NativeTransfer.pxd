#distutils: language = c++
#cython: language_level = 3
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

cimport BasicTools.Containers.NativeUnstructuredMesh as cNUM
cimport BasicTools.FE.Numberings.NativeDofNumbering as cNDN

from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool

import numpy as np
cimport numpy as cnp
cnp.import_array()

from eigency.core cimport *

from BasicTools.CythonDefs cimport CBasicIndexType, CBasicFloatType
from BasicTools.NumpyDefs import PBasicIndexType, PBasicFloatType

cdef extern from "Containers/FieldTransfer.h" namespace "BasicTools":
    cdef cppclass TransferClass:
        string ToStr()
        void SetVerbose(bool)
        void SetTransferMethod(const string& method)
        string GetTransferMethod()
        void SetSourceSpace(const string& space);
        void SetSourceMesh(cNUM.UnstructuredMesh* sourceMesh);
        void SetSourceNumbering(cNDN.DofNumbering* numbering);

        void SetTargetPoints(FlattenedMapWithOrder[Matrix, CBasicFloatType, Dynamic, Dynamic, RowMajor]& )
        void Compute() nogil
        PlainObjectBase& GetStatus();

        vector[CBasicIndexType] rows;
        vector[CBasicIndexType] cols;
        vector[CBasicFloatType] data;

        CBasicIndexType nb_source_Dofs;
        CBasicIndexType nb_targetPoints;
        bool useEdges;

cdef class NativeTransfer:

    cdef TransferClass cpp_object
    cdef TransferClass* GetCppPointer(self) nogil

    cdef cNUM.CUnstructuredMesh sourceMesh
    cdef dict __dict__

    #def SetVerbose(self, bool)
    #def SetTransferMethod(self, str)

    #def ToStr(self)
    #def ComputeShapeFunctionsOnElement(self, coordAtDofs, localspace, localnumbering, point, faceElementType)
    #def ComputeBarycentricCoordinateOnElement(self, coordAtDofs,localspace,targetPoint,elementType)
