# distutils: language = c++
#cython: language_level=3

cimport numpy as np
from libcpp.string cimport string
from libcpp.vector cimport vector

from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "FE/NativeNumericalWeakForm.h" :
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


cdef extern from "FE/NativeNumericalWeakForm.h" :
    cdef cppclass WeakMonom:
        WeakMonom() except     +
        double prefactor;
        vector[WeakTerm] prod

cdef extern from "FE/NativeNumericalWeakForm.h" :
    cdef cppclass WeakForm:
        WeakForm() except     +
        vector[WeakMonom] form
        int GetNumberOfTerms()

cdef class PyWeakTerm:
    cdef WeakTerm* _c_WeakTerm
    cdef bint pointerOwner
    cdef WeakTerm* GetCppObject(self)
    @staticmethod
    cdef PyWeakTerm create(WeakTerm* )
    cdef PyWeakTerm copy(self)


cdef class PyWeakMonom:
    cdef WeakMonom* _c_WeakMonom
    cdef bint pointerOwner
    cpdef AddProd(self,PyWeakTerm)
    cpdef int GetNumberOfProds(self)
    cdef PyWeakTerm GetProd(self, int)
    cdef WeakMonom* GetCppObject(self)
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
    cdef WeakForm* GetCppObject(self)
    @staticmethod
    cdef PyWeakForm create(WeakForm*)

cdef class PyWeakFormIter:
    cdef int index
    cdef PyWeakForm father

cdef class PyWeakMonomIter:
    cdef int index
    cdef PyWeakMonom father

