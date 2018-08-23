# distutils: language = c++


import cython
import numpy as np
cimport numpy as np
from cython.operator cimport dereference as deref

testcharacter = "'"

cimport WeakFormNumericalWrapper
from WeakFormNumericalWrapper cimport WeakTerm
from WeakFormNumericalWrapper cimport WeakMonom
from WeakFormNumericalWrapper cimport WeakForm


cdef class PyWeakTerm:

    def __cinit__(self):
       self._c_WeakTerm = new WeakFormNumericalWrapper.WeakTerm()
       self.pointerOwner = True
       if self._c_WeakTerm is NULL:
            raise MemoryError()

    def __dealloc__(self):
        if self._c_WeakTerm is not NULL and self.pointerOwner == True:
            del self._c_WeakTerm

    cdef WeakFormNumericalWrapper.WeakTerm* GetCppObject(self):
        return self._c_WeakTerm

    @staticmethod
    cdef PyWeakTerm create(WeakTerm* ptr):
        obj = <PyWeakTerm>PyWeakTerm.__new__(PyWeakTerm) # create instance without calling __init__
        obj._c_WeakTerm = ptr
        obj.pointerOwner = False
        return obj

    @property
    def derDegree(self):
        return self._c_WeakTerm.derDegree

    @derDegree.setter
    def derDegree(self, value):
        self._c_WeakTerm.derDegree = value

    @property
    def fieldName(self):
        return self._c_WeakTerm.fieldName

    @fieldName.setter
    def fieldName(self, value):
        self._c_WeakTerm.fieldName = value

    @property
    def derCoordName(self):
        return self._c_WeakTerm.derCoordName

    @derCoordName.setter
    def derCoordName(self, value):
        self._c_WeakTerm.derCoordName = value

    @property
    def derCoordIndex_(self):
        return self._c_WeakTerm.derCoordIndex_

    @derCoordIndex_.setter
    def derCoordIndex_(self, value):
        self._c_WeakTerm.derCoordIndex_ = value

    @property
    def constant(self):
        return self._c_WeakTerm.constant

    @constant.setter
    def constant(self, value):
        self._c_WeakTerm.constant = value

    @property
    def spaceIndex_(self):
        return self._c_WeakTerm.spaceIndex_

    @spaceIndex_.setter
    def spaceIndex_(self, value):
        self._c_WeakTerm.spaceIndex_ = value

    @property
    def numberingIndex_(self):
        return self._c_WeakTerm.numberingIndex_

    @numberingIndex_.setter
    def numberingIndex_(self, value):
        self._c_WeakTerm.numberingIndex_ = value

    @property
    def normal(self):
        return self._c_WeakTerm.normal

    @normal.setter
    def normal(self, value):
        self._c_WeakTerm.normal = value

    @property
    def valuesIndex_(self):
        return self._c_WeakTerm.valuesIndex_

    @valuesIndex_.setter
    def valuesIndex_(self, value):
        self._c_WeakTerm.valuesIndex_ = value


    @property
    def numberingIndex_(self):
        return self._c_WeakTerm.numberingIndex_

    @numberingIndex_.setter
    def numberingIndex_(self, value):
        self._c_WeakTerm.numberingIndex_ = value

    @property
    def internalType(self):
        return self._c_WeakTerm.internalType

    @internalType.setter
    def internalType(self, value):
        self._c_WeakTerm.internalType = value

    def __str__(self):
        res = ""
        if self.derDegree > 0 and self.normal == 0 :
#            #res += "d" + self.fieldName + "/"  + "d"  + str(self.derCoordName)
            res += "Derivative("+str(self.fieldName)+"(x, y, z), "+str(self.derCoordName)+")"
        else:
            res += self.fieldName
        return res

cdef class PyWeakMonom:
    def __cinit__(self):
        self.pointerOwner = True
        self._c_WeakMonom = new WeakFormNumericalWrapper.WeakMonom()
        if self._c_WeakMonom is NULL:
            raise MemoryError()

    def __dealloc__(self):
        if self._c_WeakMonom is not NULL and self.pointerOwner == True:
            del self._c_WeakMonom

    cpdef  hasVariable(self,str var):
        for m in self :
            if m.fieldName == str(var):
                return True
        return False

    cpdef AddProd(self,PyWeakTerm term):
        term.pointerOwner = False
        self._c_WeakMonom.prod.push_back(deref(term.GetCppObject()))

    cpdef int GetNumberOfProds(self):
        return self._c_WeakMonom.prod.size()

    cdef PyWeakTerm GetProd(self, int n):
        cdef WeakTerm* res = &self._c_WeakMonom.prod[n]
        return PyWeakTerm.create(res)

    cdef WeakFormNumericalWrapper.WeakMonom* GetCppObject(self):
        return self._c_WeakMonom

    @staticmethod
    cdef PyWeakMonom create(WeakMonom* ptr):
        obj = <PyWeakMonom>PyWeakMonom.__new__(PyWeakMonom) # create instance without calling __init__
        obj._c_WeakMonom = ptr
        obj.pointerOwner = False
        return obj

    @property
    def prefactor(self):
        return self._c_WeakMonom.prefactor

    @prefactor.setter
    def prefactor(self, value):
        self._c_WeakMonom.prefactor = value


    def __str__(self):
        res = str(self._c_WeakMonom.prefactor)
        for i in range(self.GetNumberOfProds()):
            res += "*"
            res += str(self.GetProd(i))
        return res

    def __iter__(self):
        return PyWeakMonomIter(self)

cdef class PyWeakForm:
    def __cinit__(self):
        self.pointerOwner = False
        self._c_WeakForm = new WeakForm()
        if self._c_WeakForm is NULL:
            raise MemoryError()

    def __dealloc__(self):
        if self._c_WeakForm is not NULL and self.pointerOwner == True:
            del self._c_WeakForm

    cdef WeakForm* GetCppObject(self):
        return self._c_WeakForm

    cpdef AddTerm(self,PyWeakMonom monom):
        monom.pointerOwner = False
        self._c_WeakForm.form.push_back(deref(monom.GetCppObject()))

    cpdef int GetNumberOfTerms(self):
        return self._c_WeakForm.GetNumberOfTerms()

    cpdef PyWeakMonom GetMonom(self, int n):
        cdef WeakMonom* res = &self._c_WeakForm.form[n]
        return PyWeakMonom.create(res)

    def GetRightPart(self,list unknownvars):
        res = PyWeakForm()
#        res.pointerOwner
#        #res.OnTag = self.OnTag
#        cdef WeakMonom p ;
        for p in self:
            for uv in unknownvars:
                if p.hasVariable(uv):
                    break
            else:
                res.AddTerm(p)
        return res

    def  GetLeftPart(self,list unknownvars):
        res = PyWeakForm()
#        #res.OnTag = self.OnTag
        for p in self:
            tocopy = False
            for uv in unknownvars:
               if p.hasVariable(uv):
                    tocopy =True
                    break
            if tocopy:
                res.AddTerm(p)
        return res

    def __str__(self):
        res = ""
        for i in range(self.GetNumberOfTerms()):
            res +=str(self.GetMonom(i))   + "\n"
        return res


    @staticmethod
    cdef PyWeakForm create(WeakForm* ptr):
        obj = <PyWeakForm>PyWeakForm.__new__(PyWeakForm) # create instance without calling __init__
        obj._c_WeakForm = ptr
        return obj

    def __iter__(self):
        return PyWeakFormIter(self)

cdef class PyWeakFormIter():
    def __init__(self,PyWeakForm father):
        self.index = father.GetNumberOfTerms()
        self.father = father

    def __next__(self):
        if self.index == 0:
            raise StopIteration
        self.index = self.index - 1
        return self.father.GetMonom(self.index)

    def next(self):
        return self.__next__()


cdef class PyWeakMonomIter():
    def __init__(self,PyWeakMonom father):
        self.index = father.GetNumberOfProds()
        self.father = father

    def __next__(self):
        if self.index == 0:
            raise StopIteration
        self.index = self.index - 1
        return self.father.GetProd(self.index)

    def next(self):
        return self.__next__()

def CheckIntegrity(GUI=False):
    F = PyWeakForm()
    M = PyWeakMonom()
    T = PyWeakTerm()
    M.AddProd(T)
    F.AddTerm(M)

    for term in F:
        print(term)

    return 'OK'

#@cython.final
#cdef class Weakterm(object):
#    #cdef CWeakterm *ptr
#    #cdef bint allocated
#    def __cinit__(self,CWeakterm ptr ):
#        if ptr is None:
#            self.ptr = new CWeakterm()
#            self.allocated = True
#        else:
#            self.ptr = prt
#            self.allocated = False
#
#    def __dealloc__(self):
#        if self.allocated:
#            del self.ptr
#cdef class Weakterm:
##    cdef public char* fieldName
##    cdef public char* derCoordName
##    cdef public int derDegree
##    cdef public bint constant
##    cdef public bint normal
#
##    cdef public int spaceIndex_
##    cdef public int derCoordIndex_
##    cdef public int numberingIndex_
##    cdef public int valuesIndex_
##    cdef public int internalType
#
#    def __init__(self):
#        self.fieldName = ""
#        self.derCoordName = ""
#        self.derDegree = -1
#        self.constant = False
#        self.normal = False
#
#        self.spaceIndex_ = 0
#        self.derCoordIndex_ = 0
#        self.numberingIndex_ = 0
#        self.valuesIndex_ = 0
#        self.internalType = 0
#        # 0 normal
#        # 1 constant
#        # 2 unkown
#        # 3 test
#        # 4 extra
#
#        # 5 Funciton ([x,y,x]) experimental
#
#

#
#
#@cython.final
#cdef public class WeakMonom[object CWeakMonom, type CWeakMonomType ]:
##    cdef DTYPEfloat_t prefactor
##    cdef public list prod
#
#    @property
#    def prefactor(self):
#        return self.prefactor
#
#    @prefactor.setter
#    def prefactor(self, value):
#        self.prefactor = value
#
#
#    def __init__(self):
#        self.prefactor = 1.
#        self.prod = []
#
#
#    cpdef bint hasTestFunc(self):
#        for p in self.prod :
#            if p.fieldName[-1] == testcharacter:
#                return True
#        return False
#
#    cpdef bint hasVariable(self,str var):
#        for m in self.prod :
#            if m.fieldName == str(var):
#                return True
#
#        return False
#
#
#@cython.final
#@cython.nonecheck(False)  # turn off negative index wrapping for entire function
#cdef public  class WeakForm[object CWeakForm, type CWeakFormType ]:
##    cdef public str OnTag
##    cdef public list form
#
#    def __init__(self):
#        #self.OnTag = ""
#        self.form = []
#
#
#
#    def __str__(self):
#        res = "Weak Form:\n"
#        res += "  Nb Terms : " + str(len(self.form)) +"\n"
#        fields = set()
#        const = set()
#        coords = set()
#        for p in self.form:
#            for f in p.prod:
#                if f.constant  :
#                    const.add(f.fieldName)
#                    continue
#                fields.add(f.fieldName)
#                if f.derDegree > 0:
#                    coords.add(f.derCoordName)
#        res += "  Constants : " + " ".join(const) + "\n"
#        res += "  Fields : " + " ".join(fields) + "\n"
#        res += "  Coordinates used : " + " ".join(coords) + "\n"
#        for p in self.form:
#            res += str(p) + "\n"
#        return res