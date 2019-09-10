# distutils: language = c++

cimport numpy as np
from libcpp.string cimport string
from libcpp.vector cimport vector

cimport WeakFormNumericalWrapper
from WeakFormNumericalWrapper cimport WeakTerm
from WeakFormNumericalWrapper cimport WeakMonom
from WeakFormNumericalWrapper cimport WeakForm

#cdef extern from "string" namespace "std":
#    cdef cppclass string:
#        char* c_str()
#        string(char *)

#ctypedef np.float_t DTYPEfloat_t

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

#cdef class Weakterm(object):
#    cdef CWeakterm *ptr
#    cdef bint allocated
#    cdef public char* fieldName
#    cdef public char* derCoordName
#    cdef public int derDegree
#    cdef public bint constant
#    cdef public bint normal
#
#    cdef public int spaceIndex_
#    cdef public int derCoordIndex_
#    cdef public int numberingIndex_
#    cdef public int valuesIndex_
#    cdef public int internalType




#
#cdef public class WeakMonom[object CWeakMonom, type CWeakMonomType ]:
#    cdef float prefactor
#    cdef public list prod
#    cpdef bint hasTestFunc(self)
#    cpdef bint hasVariable(self,str)
#
#
#
#cdef public class WeakForm[object CWeakForm, type CWeakFormType ]:
#    #cdef public str OnTag
#    cdef public list form
#    cpdef WeakForm GetRightPart(self,list)
#
