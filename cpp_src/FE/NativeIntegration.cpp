//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//

#include "NativeIntegration.h"
#include <Eigen/SparseCore>
#include <LinAlg/EigenTypes.h>

//http://cython.readthedocs.io/en/latest/src/userguide/wrapping_CPlusPlus.html

#include <iostream>
#include <string>

namespace BasicTools
{

//typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double

const int EnumError = -1;
const int EnumNormal = 0;
const int EnumConstant = 1;
const int EnumUnknownField = 2;
const int EnumTestField = 3;
const int EnumExtraField = 4;
const int EnumExtraIPField = 5;

MonoElementsIntegralCpp::MonoElementsIntegralCpp (): nodes(0),connectivity(0),lnumbering() {
  this->totalTestDofs = 0;
  this->totalUnkownDofs = 0;
  //this->numberOfVIJ = 0;
  this->hasnormal = false;
  this->totalvijcpt = 0;
  this->onlyUpper = false;

}
//////////////////////////////////////////
void  MonoElementsIntegralCpp::SetNumberOfValues(int i){
    for(unsigned int i=0; i < this->values.size() ; ++i){
      delete this->values[i];
    }
    this->values.resize(i,0);
};
//
void  MonoElementsIntegralCpp::SetNumberOfIPValues(int i){
    for(unsigned int i=0; i < this->ipvalues.size() ; ++i){
      delete this->ipvalues[i];
    }
    this->ipvalues.resize(i,0);
};

//
void  MonoElementsIntegralCpp::SetValueI(int i, int n, int m, double* dp){
    if(this->values[i]) delete this->values[i];
    this->values[i] = new MapMatrixDD1(dp,n,m);
}
//
void  MonoElementsIntegralCpp::SetIPValueI(int i, int n, int m, double* dp){
    if(this->ipvalues[i]) delete this->ipvalues[i];
    this->ipvalues[i] = new MapMatrixDDD(dp,n,m);
}

////////////////////  Unkown Fields  ///////////////////////////////////
void MonoElementsIntegralCpp::SetNumberOfUnkownFields(const int& n){
     this->unkownDofsOffset.resize(n,1);
     this->localUnkownDofsOffset.resize(n,1);
     this->unkownDofsNumbering.resize(n,1);
}
//
void MonoElementsIntegralCpp::SetUnkownOffset(const int& n, const int& s){
     this->unkownDofsOffset(n,0) = s ;
     //std::cout << "unkownDofsOffset : " << this->unkownDofsOffset << std::endl;
}
//
void MonoElementsIntegralCpp::SetTotalUnkownDofs(const int& n){
     this->totalUnkownDofs = n ;
}
//////////////  Test Fields  ////////////////////////////////////////////
void MonoElementsIntegralCpp::SetNumberOfTestFields(const int& n){
     this->testDofsOffset.resize(n,1);
     this->localTestDofsOffset.resize(n,1);
     this->testDofsNumbering.resize(n,1);
}
//
void MonoElementsIntegralCpp::SetTotalTestDofs(const int& n){
     this->totalTestDofs = n ;
}
//
void MonoElementsIntegralCpp::SetTestOffset(const int& n, const int& s){
     this->testDofsOffset(n,0) = s ;
}
///////////////// constants ///////////////////////////////////
void MonoElementsIntegralCpp::SetNumberOfConstants(const CBasicIndexType& n){
     this->constants.resize(n,1);
}
//
void MonoElementsIntegralCpp::SetConstants(const int& n,const double& val){
     this->constants(n,0) = val ;
}
//////////////////Prepare Fast Integration related functions //////////////////////////////

void MonoElementsIntegralCpp::AllocateWorkingElementVIJ(int size){
    //this->ev.resize(size);
    //this->ei.resize(size);
    //this->ej.resize(size);
}
//
void MonoElementsIntegralCpp::SetComputeNormal(const bool& val){
    this->hasnormal = val;
}
//
void MonoElementsIntegralCpp::SetNumberOfIntegrationPoints(const int& n){
     this->ip.resize(n,3);
     this->iw.resize(n,1);
}
//
void MonoElementsIntegralCpp::SetIntegrationPointI(const int& n,const double& w,const double& p0,const double& p1,const double& p2){
     this->iw(n,0) = w;
     this->ip(n,0) = p0;
     this->ip(n,1) = p1;
     this->ip(n,2) = p2;
}
//
void MonoElementsIntegralCpp::SetPoints(double* pd, const int& rows, const int& columns){
    if(this->nodes) delete this->nodes;
    this->nodes = new MapMatrixDDD(pd,rows,columns);
}
//
void MonoElementsIntegralCpp::SetConnectivity( CBasicIndexType* pd, const int& rows, const int& columns){
    if(this->connectivity) delete this->connectivity;
    this->connectivity = new MapMatrixIDD(pd,rows,columns);
}
//
void MonoElementsIntegralCpp::ProcessWeakForm(WeakForm* wform){
    std::vector<std::vector<WeakForm> > sortedWeakForm;
    assert(0);
}
//
template<typename T>
MatrixDDD solve(T& Jinv,MapMatrixDDD &valdphidxi){
   return Jinv.solve(valdphidxi);
}

template<>
MatrixDDD solve(MatrixDDD& Jinv,MapMatrixDDD &valdphidxi){
   return Jinv*valdphidxi;
}

void inv33(MatrixDDD& m,MatrixDDD &minv,double& det ){
        //std::cout << "m   ->   " << m<< std::endl;
        det = m(0, 0) * (m(1, 1) * m(2, 2) - m(2, 1) * m(1, 2)) -
               m(0, 1) * (m(1, 0) * m(2, 2) - m(1, 2) * m(2, 0)) +
               m(0, 2) * (m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0));

        //std::cout << "det   ->   " << det  << std::endl;
        const double invdet = 1. / det;
        //std::cout << "invdet   ->   " << invdet  << std::endl;
        //std::cout << "minv   ->   " << minv  << std::endl;
        minv(0, 0) = (m(1, 1) * m(2, 2) - m(2, 1) * m(1, 2)) * invdet;
        minv(0, 1) = (m(0, 2) * m(2, 1) - m(0, 1) * m(2, 2)) * invdet;
        minv(0, 2) = (m(0, 1) * m(1, 2) - m(0, 2) * m(1, 1)) * invdet;
        minv(1, 0) = (m(1, 2) * m(2, 0) - m(1, 0) * m(2, 2)) * invdet;
        minv(1, 1) = (m(0, 0) * m(2, 2) - m(0, 2) * m(2, 0)) * invdet;
        minv(1, 2) = (m(1, 0) * m(0, 2) - m(0, 0) * m(1, 2)) * invdet;
        minv(2, 0) = (m(1, 0) * m(2, 1) - m(2, 0) * m(1, 1)) * invdet;
        minv(2, 1) = (m(2, 0) * m(0, 1) - m(0, 0) * m(2, 1)) * invdet;
        minv(2, 2) = (m(0, 0) * m(1, 1) - m(1, 0) * m(0, 1)) * invdet;

};


void GetJackAndDet(MapMatrixDDD& valdphidxi,
                   MatrixDDD&  xcoor,
                   int& Dimensionality,
                   MatrixDDD&  Jack,
                   double&    Jdet,
                   MatrixDDD&  Jinv){

    if(Dimensionality == 0){
        Jdet = 1.;
        Jack.resize(1,1);
        Jack(0,0) = 1.;
        Jinv.resize(1,1);
        Jinv(0,0) = 1.;
        return;

    }
    Jinv.resize(Dimensionality,xcoor.cols());
    Jack = valdphidxi*xcoor;
    //std::cout << " ***** Jack ***** " << std::endl;

    switch(Dimensionality) {
      case (0): {
        Jdet = 1;
        //std::cout << "0 " << std::endl; exit(1);
        break;
      } case(1) : {
	    // we have an edge
        switch(xcoor.cols()) {
           case (1):{
              //in a 1D space
              Jdet = Jack(0, 0);
              double invdet = 1 / Jdet;
              Jinv(0,0) = invdet;
              return;
              break;
           } case (2):{
              //in a 2D space
              const double der0 = Jack(0,0); // dx
              const double der1 = Jack(0,1); // dy
              Jdet = std::sqrt(der0*der0+der1*der1);
              MatrixD21 Normal;
              Normal(0) = der1;
              Normal(1) = -der0;

              MatrixDDD m;
              m.resize(2,2),
              m(0,0) = der0;
              m(1,0) = der1;
              m.col(1) = Normal;

              Jinv.noalias() = m.inverse().transpose().col(0) ;
              return;
              break;
           } case (3):{
              const double der0 = Jack(0,0); // dx
              const double der1 = Jack(0,1); // dy
              const double der2 = Jack(0,2); // dz
              Jdet = std::sqrt(der0*der0+der1*der1+der2*der2);
              MatrixD31 N0;
              MatrixD31 N1;
              if(der0 == 0 && der1 == 0){
                 N0 <<  1., 0., 0.;
                 N1 <<  0., 1., 0.;
              } else {
                  N0 <<  der1, -der0,0.;
                  N1 <<  -der2, 0., der0;
              }
              N0 *= 1./N0.norm();
              N1 *= 1./N1.norm();

              MatrixD33 m;
              m.resize(3,3),
              m.col(0) = Jack.row(0);
              m.col(1) = N0;
              m.col(2) = N1;

              Jinv.noalias() = m.inverse().transpose().col(0) ;
              return ;
              break;
           }
        }
        Jdet = Jack.norm();

        break;
      } case(2): {
   	    // we have a surface triange quad ...

        switch(xcoor.cols()) {
          case (2):{
            // in 2D space
            Jdet = Jack(0, 0) * Jack(1, 1) - Jack(1, 0) * Jack(0, 1);
            const double invdet = 1 /Jdet;

            Jinv(0,0) =  Jack(1,1)* invdet;
            Jinv(1,0) = -Jack(1,0)* invdet;
            Jinv(0,1) = -Jack(0,1)* invdet;
            Jinv(1,1) =  Jack(0,0)* invdet;
            return;
            break;
          } case (3):{
            // in 3D space
            //compute normal
            MatrixD31 Normal = Jack.row(0).head<3>().cross(Jack.row(1).head<3>());
            Jdet = Normal.norm();
            Normal /= Jdet;


            MatrixD33 m;
            m.row(0) = Jack.row(0);
            m.row(1) = Jack.row(1);
            m.row(2) = Normal;
            Jinv.noalias() = m.inverse().block<3,2>(0,0);

            return;
            break;
          }
        }
        break;
      } case(3):{
        // we have a 3D element in 3D
        inv33(Jack,Jinv,Jdet);
        return;
      } default : {
         std::cerr << "Error in the calculation of the jacobian" << std::endl;
         exit(1);
	  }
    }


    //Jinv.compute(Jack);
    std::cerr << "Error in the calculation of the jackob" << std::endl;
    exit(0);
  };

//
void GetJackAndDet(MapMatrixDDD& valdphidxi,
                   MatrixDD3&  xcoor,
                   int& Dimensionality,
                   MatrixDDD&  Jack,
                   double&    Jdet,
                   Eigen::ColPivHouseholderQR<MatrixDDD>&  Jinv){

    Jack = valdphidxi*xcoor;


    if(Dimensionality == xcoor.cols()) {
	  Jdet = Jack.determinant();

    } else {
        switch(Dimensionality) {
        case (0): {
            Jdet = 1;
            break;
        } case(1) : {
	    // we have and edge in 2D or 3D
            Jdet = Jack.norm();
            break;
        } case(2): {
            // TODO check
            Jdet = Jack.row(0).head<3>().cross(Jack.row(1).head<3>()).norm();
            break;
        } default : {
            std::cerr << "Error in the calculation of " << std::endl;
            assert(0);
	    }
        }
    }

    Jinv.compute(Jack);

  };

void GetNormal(const int&  SpaceDim, const int& elementDim,MatrixDDD& Jack,MatrixD31& Normal){

    if( elementDim == 2 && SpaceDim == 3){
        Normal = Jack.row(0).head<3>().cross(Jack.row(1).head<3>());
        Normal /= Normal.norm();
    } else if( elementDim == 1 && SpaceDim == 2 ) {
       Normal[0] = Jack(1,0)-Jack(0,0);
       Normal[1] = Jack(1,1)-Jack(0,1);
       Normal[2] = 0;
    } else if( elementDim == 1 && SpaceDim == 3 ) {
        MatrixD31 Z1;
        Z1 <<  0, 0., 1;
        Normal = Jack.row(0).head<3>().cross(Z1);
        Normal /= Normal.norm();
    } else {
        std::cerr << "error in the normal " << std::endl;
        exit(1);
    }
};
//
void MonoElementsIntegralCpp::Integrate( WeakForm* wform, std::vector<int>& idstotreat){


    bool hasright = false;
    const CBasicIndexType NumberOfTerms = wform->GetNumberOfTerms();
    const CBasicIndexType NumberOfIntegrationPoints = static_cast<CBasicIndexType>(this->iw.rows());
    const CBasicIndexType nbcols = static_cast<CBasicIndexType>(this->connectivity->cols());
    int l1 = 0;
    int l2 = 0;
    MatrixID1 leftNumbering;
    MatrixDD1 left;
    MatrixID1 rightNumbering;
    MatrixDD1 right;

    MatrixID1 centerNumbering;
    MatrixDD1 center;
    MatrixDD1 vals;


    MatrixD31 normal;
    MatrixDDD xcoor;

    xcoor.resize(nbcols,this->nodes->cols());

    MatrixDDD ElementMatrix(maxsizelocalTestDofs,maxsizelocalUnkownDofs);

    LocalSpace& geoSpace = this->lspaces[this->geoSpaceNumber];
    const int elemdim = geoSpace.dimensionality;

    const CBasicIndexType idstotreat_s = static_cast<CBasicIndexType>(idstotreat.size());
    QRType Jinv;
    int n;
    int rightIndex=0;
    int leftIndex=0;
    for(int elem_counter =0; elem_counter< idstotreat_s; ++elem_counter){
        ElementMatrix =  MatrixDDD::Zero(maxsizelocalTestDofs,maxsizelocalUnkownDofs);

        n = idstotreat[elem_counter];
//      the coordinates of the nodes
        for( int j = 0; j < nbcols;++j){
            xcoor.row(j) = this->nodes->row((*this->connectivity)(n,j));
        }
        // for the moment we must copy the data in a for loop, the next line will work
        // with eigen 3.3.5
        //https://stackoverflow.com/questions/40074738/eigen-extracting-submatrix-from-vector-of-indices
        //auto xcoor = this->nodes(this->connectivity.row(n),Eigen::placeholders::all);

        for(int ip = 0; ip< NumberOfIntegrationPoints;  ++ip){
//            """ we recover the jacobian matrix """
            MatrixDDD Jack;
            Jack.resize(3,3);
            double Jdet;

            GetJackAndDet(*geoSpace.valdphidxi[ip],
                          xcoor,
                          geoSpace.dimensionality,
                          Jack,
                          Jdet,
                          Jinv);


            for(unsigned int s=0 ; s< this->lspaces.size(); ++s){
                this->lspaces[s].SetActiveIntegrationPoint(ip,Jinv );
            }

            if(this->hasnormal)
                GetNormal(3,elemdim,Jack,normal);

            for(int termn=0; termn<NumberOfTerms;++termn){
                WeakMonom& monom =  wform->form[termn];
                double factor = monom.prefactor;
                if(this->onlyEvaluation){
//                    # For the evaluation we only add the constribution without doing the integration
//                    # the user is responsible of dividing by the mass matrix to get the correct values
//                    # also the user can use a discontinues field to generate element surface stress (for example)
                } else {
                    factor *= Jdet;
                    factor *= iw(ip,0);
                }
                hasright = false;

                const CBasicIndexType numberOfProds = static_cast<CBasicIndexType>(monom.prod.size());

                for(CBasicIndexType prodn=0; prodn< numberOfProds; ++prodn){


                    WeakTerm& term = monom.prod[prodn];

                    if(term.internalType == EnumNormal){
                        factor *= normal(term.derDegree,0);
                        continue;
                    } else if( term.internalType == EnumConstant){
                        factor *= this->constants(term.valuesIndex_,0);
                        continue;
                    } else if( term.internalType == EnumUnknownField){
                        LocalSpace& cs = this->lspaces[term.spaceIndex_];
                        if(term.derDegree == 1){
                            right = cs.GetBxByBz().row(term.derCoordIndex_);
                        } else {
                            right = cs.GetNxNyNz();
                        }
//                        rightNumbering = this->lnumbering[term.numberingIndex_]->row(n);
//                        rightNumbering.array() += this->unkownDofsOffset(term.valuesIndex_,0);
                        hasright = true;
                        l2 = this->lspaces[term.spaceIndex_].numberOfShapeFunctions;
                        rightIndex = term.valuesIndex_;
                        continue;
                    } else if( term.internalType == EnumTestField){
                        LocalSpace& cs = this->lspaces[term.spaceIndex_];
                        if(term.derDegree == 1){
                            left = cs.GetBxByBz().row(term.derCoordIndex_);
                        } else {
                            left = cs.GetNxNyNz();
                        }
                        leftNumbering = this->lnumbering[term.numberingIndex_]->row(n);
                        leftNumbering.array() += this->testDofsOffset(term.valuesIndex_,0);
                        l1 = this->lspaces[term.spaceIndex_].numberOfShapeFunctions;
                        leftIndex = term.valuesIndex_;
                        continue;
                    } else if( term.internalType == EnumExtraField){
                        LocalSpace& cs = this->lspaces[term.spaceIndex_];
                        if(term.derDegree == 1){
                            center = cs.GetBxByBz().row(term.derCoordIndex_);
                        } else {
                            center = cs.GetNxNyNz();
                        }
                        centerNumbering = this->lnumbering[term.numberingIndex_]->row(n);
                        const CBasicIndexType ss = static_cast<CBasicIndexType>(centerNumbering.rows());
                        vals.resize(ss,1);
                        for(CBasicIndexType c=0; c < ss; ++c){
                            vals(c,0) = (*(this->values[term.valuesIndex_]))(centerNumbering(c),0);
                        }

                        factor *= center.dot(vals);
                        continue;
                    } else if( term.internalType == EnumExtraIPField){
                        if(term.derDegree == 1){
                            std::cout << "Cant compute a derivative on Integration Point Fields" << term.fieldName << std::endl;
                        }
                        factor *= (*(this->ipvalues[term.valuesIndex_]))(n,ip);
                        continue;
                    } else {
                        std::cout << "Cant treat term " << term.fieldName << std::endl;
                        exit(1);
                    }
                }
                if(factor == 0) continue;



                if(hasright){

//                    const int l = l1*l2;
//                    int l2cpt = fillcpt;
//                    for(int i=0;i<l1;++i){
//                        for(int j=0;j<l2;++j){
//                            ev[l2cpt] =  left[i]*right[j]*factor;
//                            ej[l2cpt] = rightNumbering[j];
//                            l2cpt +=1;
//                        }
//                    }
//                    l2cpt = fillcpt;
//                    for(int j=0;j<l2;++j){
//                        for(int i=0;i<l1;++i){
//                            ei[l2cpt] = leftNumbering[j];
//                            l2cpt += 1;
//                        }
//                    }
//                    fillcpt += l;
                    const CBasicIndexType loff =   this->localTestDofsOffset[leftIndex];
                    const CBasicIndexType roff =   this->localUnkownDofsOffset[rightIndex];
                    for(int i=0;i<l1;++i){
                        for(int j=0;j<l2;++j){
                            ElementMatrix(i+loff,j+roff) += left[i]*right[j]*factor;
                        }
                    }


                } else {
                    for(int cpt=0;cpt< l1; ++cpt){
                        this->F[leftNumbering[cpt]] += left[cpt]*factor;
                    }

                }
            }
        }
        for (int vi=0; vi<localTestDofsOffset.rows(); ++vi){
            leftNumbering= this->lnumbering[this->testDofsNumbering[vi]]->row(n);
            leftNumbering.array() += this->testDofsOffset(vi,0);
            for (int vj= 0; vj <localUnkownDofsOffset.rows(); ++vj){
                rightNumbering = this->lnumbering[this->unkownDofsNumbering[vj]]->row(n) ;
                rightNumbering.array() += this->unkownDofsOffset(vj,0);
                for (int i=0; i<leftNumbering.rows(); ++i){
                    for (int j=0; j<rightNumbering.rows(); ++j){
                        double val = ElementMatrix(localTestDofsOffset[vi]+i,localUnkownDofsOffset[vj]+j);
                        if (val != 0){
                            const CBasicIndexType ii = leftNumbering[i];
                            const CBasicIndexType jj = rightNumbering[j];

                            if(this->onlyUpper && ii > jj ){
                                continue;
                            }

                            vK[totalvijcpt] = val;
                            iK[totalvijcpt] = ii;
                            jK[totalvijcpt] = jj;
                            ++totalvijcpt;

                        }

                    }
                }
            }
        }
    }

};


} // namespace BasicTools
