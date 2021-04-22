#distutils: language = c++
#cython: language_level = 3
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
from libcpp.vector cimport vector

cimport numpy as cnp

cimport BasicTools.FE.WeakForms.NativeNumericalWeakForm as NNWF

cdef extern from "FE/NativeIntegration.h" namespace "BasicTools" :
    cdef cppclass MonoElementsIntegralCpp:
        MonoElementsIntegralCpp() except +
        void SetNumberOfUnkownFields(int)
        void SetUnkownOffset(int,int)
        void SetTotalUnkownDofs(int)

        void SetNumberOfTestFields(int)
        void SetTestOffset(int,int)
        void SetTotalTestDofs(int)

        void SetNumberOfConstants(int)
        void SetConstants(int,double)

        void AllocateWorkingElementVIJ(int size)

        void SetNumberOfIntegrationPoints(int)
        void SetIntegrationPointI(int, double,double,double,double)
        int totalUnkownDofs,totalTestDofs

        void SetComputeNormal(int)
        void SetDomainToTreat()
        void SetPoints(double*, int, int)
        void SetConnectivity(cnp.int_t*, int, int)

        void SetNumberOfSpaces(int)
        void InitSpaceS(const int& s,
                 const int& dim ,
                 const int& NumberOfShapeFunctions,
                 const int& numberOfIntegrationPoints )
        void SetSpaceSvalNI(const int& spaceNumber,
                      const int& integrationPoint,
                      double* pd)
        void SetSpaceSvaldphidxiI(const int& spaceNumber,
                      const int& integrationPoint,
                      double* pd)


        void SetNumberOfNumberings(int i)
        void SetNumberingI(int i, int n, int m, cnp.int_t* ip)
        void SetNumberOfValues(int i)
        void SetValueI(int i, int n, int m, double* dp)

        void SetNumberOfIPValues(int i)
        void SetIPValueI(int i, int n, int m, double* dp)

        void SetLocalOffsets(int, vector[double]&,vector[int]&,int, vector[double]&,vector[int]&)
        void Integrate(NNWF.WeakForm*, vector[int]) nogil
        int GetNumberOfUsedIvij()
        int AddToNumbefOfUsedIvij(int)
        double* vK
        cnp.int_t* iK
        cnp.int_t* jK
        double* F
        int totalvijcpt
        bint onlyEvaluation
        int  geoSpaceNumber

