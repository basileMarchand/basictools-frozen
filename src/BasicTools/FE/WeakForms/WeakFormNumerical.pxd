# distutils: language = c++
#cython: language_level=3

cimport numpy as np
from libcpp.string cimport string
from libcpp.vector cimport vector

cimport BasicTools.FE.WeakForms.WeakFormNumericalWrapper as WeakFormNumericalWrapper
from BasicTools.FE.WeakForms.WeakFormNumericalWrapper cimport WeakTerm
from BasicTools.FE.WeakForms.WeakFormNumericalWrapper cimport WeakMonom
from BasicTools.FE.WeakForms.WeakFormNumericalWrapper cimport WeakForm


cdef class PyWeakTerm:
    cdef WeakFormNumericalWrapper.WeakTerm* _c_WeakTerm
    cdef bint pointerOwner
    cdef WeakFormNumericalWrapper.WeakTerm* GetCppObject(self)
    @staticmethod
    cdef PyWeakTerm create(WeakTerm* )


cdef class PyWeakMonom:
    cdef WeakMonom* _c_WeakMonom
    cdef bint pointerOwner
    cpdef AddProd(self,PyWeakTerm)
    cpdef int GetNumberOfProds(self)
    cdef PyWeakTerm GetProd(self, int)
    cdef WeakMonom* GetCppObject(self)
    cpdef hasVariable(self,str)

    @staticmethod
    cdef PyWeakMonom create(WeakMonom*)


cdef class PyWeakForm:
    cdef WeakFormNumericalWrapper.WeakForm* _c_WeakForm
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

