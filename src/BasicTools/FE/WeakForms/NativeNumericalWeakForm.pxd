# distutils: language = c++
#cython: language_level=3
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
cimport numpy as np
from libcpp.string cimport string
from libcpp.vector cimport vector
from BasicTools.CythonDefs cimport CBasicIndexType

cdef extern from "FE/NativeNumericalWeakForm.h" namespace "BasicTools" :
    cdef cppclass WeakTerm:
        WeakTerm() except +
        string fieldName
        string derCoordName
        int derDegree
        bint constant
        bint normal
        int spaceIndex_
        int derCoordIndex_
        int numberingIndex_
        int valuesIndex_
        int internalType
        int spaceIndex_
        int modeIndex_

cdef extern from "FE/NativeNumericalWeakForm.h" namespace "BasicTools" :
    cdef cppclass WeakMonom:
        WeakMonom() except +
        double prefactor
        vector[WeakTerm] prod

cdef extern from "FE/NativeNumericalWeakForm.h" namespace "BasicTools" :
    cdef cppclass WeakForm:
        WeakForm() except +
        vector[WeakMonom] form
        CBasicIndexType GetNumberOfTerms()

cdef class PyWeakTerm:
    cdef WeakTerm* cpp_object_pointer
    cdef bint pointerOwner
    cdef WeakTerm* GetCppPointer(self) nogil
    cdef void SetPointer(self, WeakTerm* obj ) nogil

cdef class PyWeakMonom:
    cdef WeakMonom* cpp_object_pointer
    cdef bint pointerOwner
    cdef WeakMonom* GetCppPointer(self) nogil
    cdef void SetPointer(self, WeakMonom* obj )

    cpdef AddProd(self, PyWeakTerm)
    cpdef CBasicIndexType GetNumberOfProds(self)
    cdef PyWeakTerm GetProd(self, CBasicIndexType)

    cpdef hasVariable(self, str)

cdef class PyWeakForm:
    cdef WeakForm* cpp_object_pointer
    cdef bint pointerOwner
    cdef WeakForm* GetCppPointer(self) nogil

    cpdef CBasicIndexType GetNumberOfTerms(self)
    cpdef PyWeakMonom GetMonom(self, int n)
    cpdef AddTerm(self, PyWeakMonom)

cdef int TN_PyWeakTerm = 0
cdef int TN_PyWeakMonom = 0
cdef int TN_PyWeakForm = 0