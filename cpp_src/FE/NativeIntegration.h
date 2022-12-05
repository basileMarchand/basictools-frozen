//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//
#pragma once

#include <string>
#include <iostream>
#include <vector>

#include <FE/NativeNumericalWeakForm.h>
#include <LinAlg/EigenTypes.h>
namespace BasicTools
{

template<typename T>
MatrixDDD solve(T& Jinv,MapMatrixDDD &valdphidxi);
template<>
MatrixDDD solve(MatrixDDD& Jinv,MapMatrixDDD &valdphidxi);

struct LocalSpace{
   int dimensionality;
   int numberOfShapeFunctions;
   int numberOfIntegrationPoints;
   int activeIntegrationPoint;

   std::vector< MapMatrixDDD* > valN;
   std::vector< MapMatrixDDD* > valdphidxi;
   MatrixDDD BxByBz;

   void Init(const int& dim ,
             const int& NumberOfShapeFunctions,
             const int& numberOfIntegrationPoints){
      this->dimensionality = dim;
      this->numberOfShapeFunctions = NumberOfShapeFunctions;
      this->activeIntegrationPoint=-1;
      this->valN.resize(numberOfIntegrationPoints,0);
      this->valdphidxi.resize(numberOfIntegrationPoints,0);
   }
   void resize(int& s){
       this->RelseaseData();
       for(unsigned i=0; i < this->valN.size() ; ++i){
           this->valN.push_back(0 );
           this->valdphidxi.push_back(0 );
       }
   }
   void SetvalNI(const int& integrationPoint,double* pd){
       this->valN[integrationPoint] =   new MapMatrixDDD(pd,this->numberOfShapeFunctions,1);
   }
   void RelseaseData(){
       for(unsigned i=0; i < this->valN.size() ; i++){
        if(this->valN[i]) delete this->valN[i];
        if(this->valdphidxi[i]) delete this->valdphidxi[i];
       }
       this->valN.resize(0);
       this->valdphidxi.resize(0);
   }
   LocalSpace(){
       this->valN.resize(0);
       this->valdphidxi.resize(0);
   }
   ~LocalSpace(){
       this->RelseaseData();
   }
   void SetvaldphidxiI(const int& integrationPoint,double* pd){
    this->valdphidxi[integrationPoint] = new MapMatrixDDD(pd,this->dimensionality,this->numberOfShapeFunctions);

   }


   MapMatrixDDD& GetNxNyNz(){
       return *this->valN[this->activeIntegrationPoint];
   }

   MatrixDDD & GetBxByBz(){
           return this->BxByBz;
   }

   void SetActiveIntegrationPoint(const int& ip,QRType& Jinv ){
       this->activeIntegrationPoint = ip;
       if(this->valdphidxi[ip])
           this->BxByBz = solve(Jinv,*this->valdphidxi[ip]);
   }
};

//
struct MonoElementsIntegralCpp{
  int totalUnkownDofs;
  MatrixID1 unkownDofsOffset;
  MatrixID1 unkownDofsNumbering;

  MatrixID1 localUnkownDofsOffset;
  int maxsizelocalUnkownDofs;

  int totalTestDofs;
  MatrixID1 testDofsOffset;
  MatrixID1 testDofsNumbering;

  MatrixID1 localTestDofsOffset;
  int maxsizelocalTestDofs;

  MatrixDD1 constants;

  MatrixDD1 iw;
  MatrixDD3 ip;

  int totalvijcpt;
  CBasicFloatType *vK;
  CBasicIndexType *iK;
  CBasicIndexType *jK;
  CBasicFloatType *F;

  bool hasnormal;
  bool onlyEvaluation;
  bool onlyUpper;

  MapMatrixDDD* nodes;
  MapMatrixIDD* connectivity;

  MonoElementsIntegralCpp();
  ~MonoElementsIntegralCpp(){
      if(this->nodes) delete this->nodes;
      if(this->connectivity) delete this->connectivity;
      for(unsigned int i=0; i < this->lnumbering.size() ; ++i){
        delete this->lnumbering[i];
      }
      for(unsigned int i=0; i < this->values.size() ; ++i){
        delete this->values[i];
      }
  }
  int geoSpaceNumber;

  std::vector<LocalSpace> lspaces;
  std::vector<MapMatrixIDD* >  lnumbering;
  std::vector<MapMatrixDD1*>  values;
  std::vector<MapMatrixDDD*>  ipvalues;

  void SetLocalOffsets(const int& maxSizeUDof,
                       const std::vector<int>& ludof,
                       const std::vector<int>& luNumberingindex,
                       const int& maxSizeTDof,
                       const std::vector<int>& ltdof,
                       const std::vector<int>& ltNumberingindex){
      this->maxsizelocalUnkownDofs = maxSizeUDof;
      for(unsigned int i =0; i < ludof.size(); ++i){
          this->localUnkownDofsOffset[i] = ludof[i];
          this->unkownDofsNumbering[i] =luNumberingindex[i];
          }
      this->maxsizelocalTestDofs = maxSizeTDof;
      for(unsigned int i =0; i < ltdof.size(); ++i){
          this->localTestDofsOffset[i] =ltdof[i];
          this->testDofsNumbering[i] =ltNumberingindex[i];
          }
  }

  void Reset(){
    this->hasnormal = false;
    this->totalvijcpt = 0;
    this->onlyUpper = false;
  };

  void SetGeoSpace(const int& i){
          this->geoSpaceNumber = i ;
  };
  void SetNumberOfSpaces(const int& i){
          this->lspaces.resize(i);
  };

  void InitSpaceS(const int& s,
                 const int& dim ,
                 const int& NumberOfShapeFunctions,
                 const int& numberOfIntegrationPoints ){
    this->lspaces[s].Init(dim, NumberOfShapeFunctions,numberOfIntegrationPoints  );
  };
  void SetSpaceSvalNI(const int& spaceNumber,
                      const int& integrationPoint,
                      double* pd){
     this->lspaces[spaceNumber].SetvalNI(integrationPoint,pd);
  }
  void SetSpaceSvaldphidxiI(const int& spaceNumber,
                      const int& integrationPoint,
                      double* pd){
     this->lspaces[spaceNumber].SetvaldphidxiI(integrationPoint, pd);
  }
  //////////////////////////////////////////
  void SetNumberOfNumberings(int i){

      for(unsigned int i=0; i < this->lnumbering.size() ; ++i){
        if(this->lnumbering[i]) {
            delete this->lnumbering[i];
            this->lnumbering[i] = nullptr;
        }
      }
      this->lnumbering.resize(i,0);
      //for(unsigned int i=0; i < this->lnumbering.size() ; ++i){
      //  this->lnumbering[i] = nullptr;
      //}

  };
  //
  void SetNumberingI(int i, int n, int m, CBasicIndexType* ip){
     if(this->lnumbering[i] != nullptr) {
         delete this->lnumbering[i];
         this->lnumbering[i] = nullptr;
     }
     this->lnumbering[i] = new MapMatrixIDD(ip,n,m);
   }
  //////////////////////////////////////////
  void SetNumberOfValues(int i);
  void SetValueI(int i, int n, int m, double* dp);
  void SetNumberOfIPValues(int i);
  void SetIPValueI(int i, int n, int m, double* dp);
  ///////////////  Unkown Fields  ///////////////////////////////////////
  void SetNumberOfUnkownFields(const int& n);
  void SetUnkownOffset(const int& n, const int& s);
  void SetTotalUnkownDofs(const int& n);
  ///////////////  Test Fields  //////////////////////////////////
  void SetNumberOfTestFields(const int& n);
  void SetTotalTestDofs(const int& n);
  void SetTestOffset(const int& n, const int& s);
  ///////////////// constants ///////////////////////////////////
  void SetNumberOfConstants(const CBasicIndexType& n);
  void SetConstants(const int& n,const double& val);
  //////////////// Working Memory ///////////////////////////////
  void AllocateWorkingElementVIJ(int size);
  ///////// PrepareFastIntegration /////////////////////////////
  void SetComputeNormal(const bool& val);
  ///////////////// Integration ///////////////////////////////////
  void SetNumberOfIntegrationPoints(const int& n);
  void SetIntegrationPointI(const int& n,const double& w,const double& p0,const double& p1,const double& p2);
  void SetPoints(double* pd, const int& rows, const int& columns);
  void SetConnectivity(CBasicIndexType* pd, const int& rows, const int& columns);
  void ProcessWeakForm(WeakForm* wform);
  void Integrate( WeakForm* wform, const CBasicIndexType& size, CBasicIndexType* pidstotreat );

};

} // namespace BasicTools
