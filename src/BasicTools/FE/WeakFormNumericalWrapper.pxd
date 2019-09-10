# distutils: language = c++

from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "src_cpp/NativeIntegration.h" :
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

cdef extern from "src_cpp/NativeIntegration.h" :
    cdef cppclass WeakMonom:
        WeakMonom() except     +
        double prefactor;
        vector[WeakTerm] prod

cdef extern from "src_cpp/NativeIntegration.h" :
    cdef cppclass WeakForm:
        WeakForm() except     +
        vector[WeakMonom] form
        int GetNumberOfTerms()


