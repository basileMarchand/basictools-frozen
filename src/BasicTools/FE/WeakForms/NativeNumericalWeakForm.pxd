# distutils: language = c++
#cython: language_level=3
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
cimport numpy as np
from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "FE/NativeNumericalWeakForm.h" namespace "BasicTools" :
    cdef cppclass WeakTerm:
        WeakTerm() except     +
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
        WeakMonom() except     +
        double prefactor;
        vector[WeakTerm] prod

cdef extern from "FE/NativeNumericalWeakForm.h" namespace "BasicTools" :
    cdef cppclass WeakForm:
        WeakForm() except     +
        vector[WeakMonom] form
        int GetNumberOfTerms()

cdef class PyWeakTerm:
    cdef WeakTerm* _c_WeakTerm
    cdef bint pointerOwner
    cdef WeakTerm* GetCppPointer(self)
    @staticmethod
    cdef PyWeakTerm create(WeakTerm* )
    cdef PyWeakTerm copy(self)


cdef class PyWeakMonom:
    cdef WeakMonom* _c_WeakMonom
    cdef bint pointerOwner
    cpdef AddProd(self,PyWeakTerm)
    cpdef int GetNumberOfProds(self)
    cdef PyWeakTerm GetProd(self, int)
    cdef WeakMonom* GetCppPointer(self)
    cpdef hasVariable(self,str)
    #def PyWeakMonom copy(self)

    @staticmethod
    cdef PyWeakMonom create(WeakMonom*)


cdef class PyWeakForm:
    cdef WeakForm* _c_WeakForm
    cdef bint pointerOwner
    cpdef int GetNumberOfTerms(self)
    cpdef PyWeakMonom GetMonom(self, int n)
    cpdef AddTerm(self,PyWeakMonom)
    cdef WeakForm* GetCppPointer(self)
    @staticmethod
    cdef PyWeakForm create(WeakForm*)

cdef class PyWeakFormIter:
    cdef int index
    cdef PyWeakForm father

cdef class PyWeakMonomIter:
    cdef int index
    cdef PyWeakMonom father

