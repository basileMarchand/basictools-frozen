# distutils: language = c++
#cython: language_level=3
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
from typing import Optional
import cython
import numpy as np
cimport numpy as np
from cython.operator cimport dereference as deref

testcharacter = "'"

cimport BasicTools.FE.WeakForms.NativeNumericalWeakForm as NNWF
from BasicTools.CythonDefs cimport CBasicIndexType



cdef class PyWeakTerm:
    cdef NNWF.WeakTerm* GetCppPointer(self) nogil:
        return self.cpp_object_pointer

    def __cinit__(self):
       self.cpp_object_pointer = new NNWF.WeakTerm()
       self.pointerOwner = True
       NNWF.TN_PyWeakTerm += 1

    def __dealloc__(self):
        if self.pointerOwner == True:
            del self.cpp_object_pointer
            NNWF.TN_PyWeakTerm -= 1
        #print(f"{NNWF.TN_PyWeakTerm=}")

    cdef void SetPointer(self, NNWF.WeakTerm* obj ) nogil:
        if self.pointerOwner == True:
            del self.cpp_object_pointer
        self.cpp_object_pointer = obj
        self.pointerOwner = False

    @property
    def derDegree(self):
        return deref(self.cpp_object_pointer).derDegree

    @derDegree.setter
    def derDegree(self, value):
        deref(self.cpp_object_pointer).derDegree = value

    @property
    def fieldName(self):
        return deref(self.cpp_object_pointer).fieldName.decode()

    @fieldName.setter
    def fieldName(self, value):
        deref(self.cpp_object_pointer).fieldName = value.encode()

    @property
    def derCoordName(self):
        return deref(self.cpp_object_pointer).derCoordName.decode()

    @derCoordName.setter
    def derCoordName(self, value):
        deref(self.cpp_object_pointer).derCoordName = value.encode()

    @property
    def derCoordIndex_(self):
        return deref(self.cpp_object_pointer).derCoordIndex_

    @derCoordIndex_.setter
    def derCoordIndex_(self, value):
        deref(self.cpp_object_pointer).derCoordIndex_ = value

    @property
    def constant(self):
        return deref(self.cpp_object_pointer).constant

    @constant.setter
    def constant(self, value):
        deref(self.cpp_object_pointer).constant = value

    @property
    def spaceIndex_(self):
        return deref(self.cpp_object_pointer).spaceIndex_

    @spaceIndex_.setter
    def spaceIndex_(self, value):
        deref(self.cpp_object_pointer).spaceIndex_ = value

    @property
    def numberingIndex_(self):
        return deref(self.cpp_object_pointer).numberingIndex_

    @numberingIndex_.setter
    def numberingIndex_(self, value):
        deref(self.cpp_object_pointer).numberingIndex_ = value

    @property
    def normal(self):
        return deref(self.cpp_object_pointer).normal

    @normal.setter
    def normal(self, value):
        deref(self.cpp_object_pointer).normal = value

    @property
    def valuesIndex_(self):
        return deref(self.cpp_object_pointer).valuesIndex_

    @valuesIndex_.setter
    def valuesIndex_(self, value):
        deref(self.cpp_object_pointer).valuesIndex_ = value

    @property
    def modeIndex_(self):
        return deref(self.cpp_object_pointer).modeIndex_

    @modeIndex_.setter
    def modeIndex_(self, value):
        deref(self.cpp_object_pointer).modeIndex_ = value

    @property
    def numberingIndex_(self):
        return deref(self.cpp_object_pointer).numberingIndex_

    @numberingIndex_.setter
    def numberingIndex_(self, value):
        deref(self.cpp_object_pointer).numberingIndex_ = value

    @property
    def internalType(self):
        return deref(self.cpp_object_pointer).internalType

    @internalType.setter
    def internalType(self, value):
        deref(self.cpp_object_pointer).internalType = value

    def __str__(self):
        res = ""
        if self.derDegree > 0 and self.normal == 0 :
            res += "Derivative("+str(self.fieldName)+", "+str(self.derCoordName)+")"
        else:
            res += self.fieldName
        return res

    @property
    def EnumError(self):
        return -1

    @property
    def EnumNormal(self):
        return 0

    @property
    def EnumConstant(self):
        return 1

    @property
    def EnumUnknownField(self):
        return 2

    @property
    def EnumTestField(self):
        return 3

    @property
    def EnumExtraField(self):
        return 4

    @property
    def EnumExtraIPField(self):
        return 5

cdef class PyWeakMonom:
    cdef NNWF.WeakMonom* GetCppPointer(self) nogil:
        return self.cpp_object_pointer

    def __cinit__(self):
       self.cpp_object_pointer = new NNWF.WeakMonom()
       self.pointerOwner = True
       NNWF.TN_PyWeakMonom += 1
      #print(f"alloc {NNWF.TN_PyWeakMonom=}")

    def __dealloc__(self):
        if self.pointerOwner == True:
            deref(self.cpp_object_pointer).prod.resize(0)
            del self.cpp_object_pointer
            NNWF.TN_PyWeakMonom -= 1
        #print(f"sealloc {self.pointerOwner}  {NNWF.TN_PyWeakMonom=} ")

    cdef void SetPointer(self, NNWF.WeakMonom* obj ):
        #print("SetPointer")
        if self.pointerOwner == True:
            deref(self.cpp_object_pointer).prod.resize(0)
            del self.cpp_object_pointer

        self.cpp_object_pointer = obj
        self.pointerOwner = False

    cpdef AddProd(self, PyWeakTerm term):
        deref(self.cpp_object_pointer).prod.push_back(deref(term.GetCppPointer()))

    cpdef CBasicIndexType GetNumberOfProds(self):
        return deref(self.cpp_object_pointer).prod.size()

    cdef PyWeakTerm GetProd(self, CBasicIndexType n):
        res = PyWeakTerm()
        res.SetPointer(&(deref(self.cpp_object_pointer).prod[n]))
        return res

    cpdef hasVariable(self,str var):
        for m in self :
            if m.fieldName == str(var):
                return True
        return False

    @property
    def prefactor(self):
        return deref(self.cpp_object_pointer).prefactor

    @prefactor.setter
    def prefactor(self, value):
        deref(self.cpp_object_pointer).prefactor = value

    def __iter__(self):
        return (self.GetProd(i) for i in range(self.GetNumberOfProds()))

cdef class PyWeakForm:
    def __cinit__(self):
        self.cpp_object_pointer = new NNWF.WeakForm()
        self.pointerOwner = True
        NNWF.TN_PyWeakForm += 1
        #print(f"alloc {NNWF.TN_PyWeakForm=}")

    def __dealloc__(self):
        if self.pointerOwner == True:
            deref(self.cpp_object_pointer).form.resize(0)
            del self.cpp_object_pointer
        NNWF.TN_PyWeakForm -= 1
        #print(f"alloc {NNWF.TN_PyWeakForm=}")

    cdef NNWF.WeakForm* GetCppPointer(self) nogil:
        return self.cpp_object_pointer

    cpdef AddTerm(self, PyWeakMonom monom):
        deref(self.cpp_object_pointer).form.push_back(deref(monom.GetCppPointer()))

    cpdef CBasicIndexType GetNumberOfTerms(self):
        return deref(self.cpp_object_pointer).GetNumberOfTerms()

    cpdef PyWeakMonom GetMonom(self, int n):
        res = PyWeakMonom()
        res.SetPointer(&(deref(self.cpp_object_pointer).form[n]))
        return res

    def GetRightPart(self,list unknownvars):
        res = PyWeakForm()
        for p in self:
            for uv in unknownvars:
                if p.hasVariable(uv):
                    break
            else:
                res.AddTerm(p)
        return res

    def  GetLeftPart(self,list unknownvars):
        res = PyWeakForm()
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

    def __iter__(self):
        return (self.GetMonom(i) for i in range(self.GetNumberOfTerms()))
