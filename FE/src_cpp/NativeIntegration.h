#ifndef NATIVEINTEGRATION_H
#define NATIVEINTEGRATION_H
#include <string>
#include <iostream>
#include <vector>

#include <FE/src_cpp/WeakFormNumerical.h>
#include <FE/src_cpp/Eigentype.h>

#define IN_DEBUG
#ifdef IN_DEBUG
  #define PRINTL(var) std::cout  << __FILE__ <<"::" << __LINE__ << " " << #var  << " -> " << (var) << std::endl;
  #define PRINT(var) std::cout  << #var  << " -> " << (var) << std::endl;
  #define PRINTS(str) std::cout <<  #str   << std::endl;
  template<typename T>
  void printMatrix(T& mat){
    std::cout << std::endl <<" c++  ";
    for(int i =0; i < mat.rows(); ++i){
        std::cout << std::endl;
        for(int j =0; j < mat.cols(); ++j){
            std::cout <<  mat(i,j) << " ";
        };
    };
    std::cout << std::endl;
  };
#else
  template<typename T>
  void printMatrix(T& mat){};
  #define PRINTL(var)
  #define PRINT(var)
  #define PRINTS(str)
#endif




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
       //new (this->valN[integrationPoint]) Eigen::Map<MatrixDD1>(pd,this->numberOfShapeFunctions,1);
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
       //PRINT(integrationPoint);
       //PRINT(this->dimensionality);
       //PRINT(this->numberOfShapeFunctions);
       //new (this->valdphidxi[integrationPoint]) Eigen::Map<MatrixDDD>(pd,this->dimensionality,this->numberOfShapeFunctions);
       //PRINT(this->valdphidxi[integrationPoint]);

       //printMatrix(*this->valdphidxi[integrationPoint]);
   }


   MapMatrixDDD& GetNxNyNz(){
       return *this->valN[this->activeIntegrationPoint];
   }

   MatrixDDD & GetBxByBz(){
           return this->BxByBz;
   }

   //void SetNumberOfIntegrationPoints(const int& n ){
   //     this->valN.resize(n);
   //     this->valdphidxi.resize(n);
   //}

   void SetActiveIntegrationPoint(const int& ip,QRType& Jinv ){
       this->activeIntegrationPoint = ip;
       /*
       cvcStorage2D<double> res;
        res.resize(GetNPhi(),IntegrationWeights.size() );

        for(int p = 0; p < IntegrationWeights.GetSize(); ++p) {
            //Eigen::ColPivHouseholderQR<cvcMatrixD> Jinv = jacobian(p).colPivHouseholderQr();
            for(int i = 0; i < dPhidXiAtIntegrationPoints(p).cols(); ++i) {
                Eigen::Matrix<double,Eigen::Dynamic,1> tmp  =  Jinv(p).solve(dPhidXiAtIntegrationPoints(p).col(i)) ;
                res(i,p) = tmp(derivation);
            }
        }
        return res;
        */
//       std::cout << "this->valdphidxi[ip] : " << std::endl  << *this->valdphidxi[ip] << std::endl;;


       this->BxByBz = Jinv.solve(*this->valdphidxi[ip]);


//       std::cout << "BxByBz["<< ip <<  "]"  << std::endl << this->BxByBz<< std::endl ;
   }

};
//

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
  FLOAT_TYPE *vK;
  INT_TYPE *iK;
  INT_TYPE *jK;
  FLOAT_TYPE *F;

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
  //LocalSpace lgeoSpace;
  int geoSpaceNumber;

  std::vector<LocalSpace> lspaces;
  std::vector<MapMatrixIDD* >  lnumbering;
  std::vector<MapMatrixDD1*>  values;

  //std::vector<FLOAT_TYPE> ev;
  //std::vector<INT_TYPE> ei;
  //std::vector<INT_TYPE> ej;

  void SetLocalOffsets(const int& maxSizeUDof,
                       const std::vector<double>& ludof,
                       const std::vector<int>& luNumberingindex,
                       const int& maxSizeTDof,
                       const std::vector<double>& ltdof,
                       const std::vector<int>& ltNumberingindex){
      this->maxsizelocalUnkownDofs = maxSizeUDof;
      for(unsigned int i =0; i < ludof.size(); ++i){
          this->localUnkownDofsOffset[i] =ludof[i];
          this->unkownDofsNumbering[i] =luNumberingindex[i];
          }
      this->maxsizelocalTestDofs = maxSizeTDof;
      for(unsigned int i =0; i < ltdof.size(); ++i){
          this->localTestDofsOffset[i] =ltdof[i];
          this->testDofsNumbering[i] =ltNumberingindex[i];
          }


  }

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
    //std::cout << "init Space " << s << " "
    //                           << dim <<" "
    //                           << NumberOfShapeFunctions << " "
    //                           << numberOfIntegrationPoints <<std::endl;
    this->lspaces[s].Init(dim, NumberOfShapeFunctions,numberOfIntegrationPoints  );
  };
  void SetSpaceSvalNI(const int& spaceNumber,
                      const int& integrationPoint,
                      double* pd){
     //std::cout << "SetSpaceSvalNI " << spaceNumber << " " << integrationPoint << std::endl;
     this->lspaces[spaceNumber].SetvalNI(integrationPoint,pd);
     //new (&this->lspaces[spaceNumber].valN[integrationPoint]) Eigen::Map<MatrixDD1>(pd,this->lspaces[spaceNumber].numberOfShapeFunctions,1);
  }
  void SetSpaceSvaldphidxiI(const int& spaceNumber,
                      const int& integrationPoint,
                      double* pd){
     this->lspaces[spaceNumber].SetvaldphidxiI(integrationPoint, pd);
     //new (&this->lspaces[spaceNumber].valdphidxi[integrationPoint]) Eigen::Map<MatrixDDD>(pd,this->lspaces[spaceNumber].dimensionality,this->lspaces[spaceNumber].numberOfShapeFunctions);
  }
  //////////////////////////////////////////
  void SetNumberOfNumberings(int i){
      for(unsigned int i=0; i < this->lnumbering.size() ; ++i){
        if(this->lnumbering[i])
            delete this->lnumbering[i];
      }
      this->lnumbering.resize(i,0);
  };
  //
  void SetNumberingI(int i, int n, int m, INT_TYPE* ip){
     if(this->lnumbering[i]) delete this->lnumbering[i];
     this->lnumbering[i] = new MapMatrixIDD(ip,n,m);
     //printMatrix(*this->lnumbering[i]);
   }
  //////////////////////////////////////////
  void SetNumberOfValues(int i){
      for(unsigned int i=0; i < this->values.size() ; ++i){
        delete this->values[i];
      }
      this->values.resize(i,0);
  };
  //
  void SetValueI(int i, int n, int m, double* dp){
      if(this->values[i]) delete this->values[i];
      this->values[i] = new MapMatrixDD1(dp,n,m);
   }

  ///////////////  Unkown Fields  ///////////////////////////////////////
  void SetNumberOfUnkownFields(const int& n);
  void SetUnkownOffset(const int& n, const int& s);
  void SetTotalUnkownDofs(const int& n);
  ///////////////  Test Fields  //////////////////////////////////
  void SetNumberOfTestFields(const int& n);
  void SetTotalTestDofs(const int& n);
  void SetTestOffset(const int& n, const int& s);
  ///////////////// constants ///////////////////////////////////
  void SetNumberOfConstants(const int& n);
  void SetConstants(const int& n,const double& val);
  //////////////// Working Memory ///////////////////////////////
  void AllocateWorkingElementVIJ(int size);
  ///////// PrepareFastIntegration /////////////////////////////
    void SetComputeNormal(const bool& val);
  ///////////////// Integration ///////////////////////////////////
  void SetNumberOfIntegrationPoints(const int& n){
         this->ip.resize(n,3);
         this->iw.resize(n,1);
     //    this->geoSpace.SetNumberOfIntegrationPoints(n);
     //#    for(int i =0, i < this->__usedSpaces__.size(); ++i){
     //       this->__usedSpaces__.SetNumberOfIntegrationPoints(n);
     //   }
  }
  void SetIntegrationPointI(const int& n,const double& w,const double& p0,const double& p1,const double& p2){
      this->iw(n,0) = w;
      this->ip(n,0) = p0;
      this->ip(n,1) = p1;
      this->ip(n,2) = p2;
  }
  void SetPoints(double* pd, const int& rows, const int& columns);
  //
  void SetConnectivity(INT_TYPE* pd, const int& rows, const int& columns);

  void ProcessWeakForm(WeakForm* wform){
    std::vector<std::vector<WeakForm> > sortedWeakForm;
    assert(0);
  }
  void Integrate( WeakForm* wform, std::vector<int>& idstotreat);

};

#endif