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
        void Reset()
        void SetNumberOfUnkownFields(CBasicIndexType)
        void SetUnkownOffset(int,int)
        void SetTotalUnkownDofs(int)

        void SetNumberOfTestFields(CBasicIndexType)
        void SetTestOffset(int,int)
        void SetTotalTestDofs(int)

        void SetNumberOfConstants(CBasicIndexType)
        void SetConstants(int,double)

        void AllocateWorkingElementVIJ(int size)

        void SetNumberOfIntegrationPoints(int)
        void SetIntegrationPointI(int, double,double,double,double)
        int totalUnkownDofs,totalTestDofs

        void SetComputeNormal(int)
        void SetDomainToTreat()
        void SetPoints(double*, int, int)
        void SetConnectivity(CBasicIndexType*, int, int)

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
        void SetNumberingI(int i, int n, int m, CBasicIndexType* ip)
        void SetNumberOfValues(CBasicIndexType i)
        void SetValueI(int i, int n, int m, double* dp)

        void SetNumberOfIPValues(int i)
        void SetIPValueI(int i, int n, int m, double* dp)

        void SetLocalOffsets(int, vector[int]&,vector[int]&,int, vector[int]&,vector[int]&)
        void Integrate(NNWF.WeakForm*, vector[int]) nogil
        int GetNumberOfUsedIvij()
        int AddToNumbefOfUsedIvij(int)
        double* vK
        CBasicIndexType* iK
        CBasicIndexType* jK
        double* F
        int totalvijcpt
        bint onlyEvaluation
        int  geoSpaceNumber

