#distutils: language = c++
#cython: language_level = 3
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
from libcpp.vector cimport vector
from BasicTools.CythonDefs cimport CBasicIndexType

cimport BasicTools.FE.WeakForms.NativeNumericalWeakForm as NNWF

cdef extern from "FE/NativeIntegration.h" namespace "BasicTools" :
    cdef cppclass MonoElementsIntegralCpp:
        MonoElementsIntegralCpp() except +
        void Reset() nogil
        void SetNumberOfUnkownFields(CBasicIndexType) nogil
        void SetUnkownOffset(int,int)  nogil
        void SetTotalUnkownDofs(int)  nogil

        void SetNumberOfTestFields(CBasicIndexType) nogil
        void SetTestOffset(int,int) nogil
        void SetTotalTestDofs(int) nogil

        void SetNumberOfConstants(CBasicIndexType) nogil
        void SetConstants(int,double) nogil

        void AllocateWorkingElementVIJ(int size) nogil

        void SetNumberOfIntegrationPoints(int) nogil
        void SetIntegrationPointI(int, double,double,double,double) nogil
        int totalUnkownDofs,totalTestDofs

        void SetComputeNormal(int) nogil
        void SetDomainToTreat() nogil
        void SetPoints(double*, int, int) nogil
        void SetConnectivity(CBasicIndexType*, int, int) nogil

        void SetNumberOfSpaces(int) nogil
        void InitSpaceS(const int& s,
                 const int& dim ,
                 const int& NumberOfShapeFunctions,
                 const int& numberOfIntegrationPoints ) nogil
        void SetSpaceSvalNI(const int& spaceNumber,
                      const int& integrationPoint,
                      double* pd) nogil
        void SetSpaceSvaldphidxiI(const int& spaceNumber,
                      const int& integrationPoint,
                      double* pd) nogil


        void SetNumberOfNumberings(int i) nogil
        void SetNumberingI(int i, int n, int m, CBasicIndexType* ip) nogil
        void SetNumberOfValues(CBasicIndexType i) nogil
        void SetValueI(int i, int n, int m, double* dp) nogil

        void SetNumberOfIPValues(int i) nogil
        void SetIPValueI(int i, int n, int m, double* dp) nogil

        void SetLocalOffsets(int, vector[int]&, vector[int]&, int, vector[int]&, vector[int]&) nogil
        void Integrate(NNWF.WeakForm*, CBasicIndexType idsize, CBasicIndexType* ids )  nogil
        int GetNumberOfUsedIvij() nogil
        int AddToNumbefOfUsedIvij(int)  nogil
        double* vK
        CBasicIndexType* iK
        CBasicIndexType* jK
        double* F
        int totalvijcpt
        bint onlyEvaluation
        int  geoSpaceNumber

