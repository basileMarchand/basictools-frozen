//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//

#include "NativeIntegration.h"
#include "LinAlg/BasicOperations.h"
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
    ElementMatrix =  MatrixDDD::Zero(maxsizelocalTestDofs,maxsizelocalUnkownDofs);
    MatrixDDD Jack(3,3);

    for(int elem_counter =0; elem_counter< idstotreat_s; ++elem_counter){
        ElementMatrix.setZero();
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

            double Jdet;

            GetInv_Jacobian_Det(*geoSpace.valdphidxi[ip],
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
                if(factor == 0.0) continue;

                if(hasright){

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
                        const double val = ElementMatrix(localTestDofsOffset[vi]+i,localUnkownDofsOffset[vj]+j);
                        if (val != 0.0){
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
