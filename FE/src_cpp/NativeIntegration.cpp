#include <iostream>

#include "NativeIntegration.h"
#include <SparseCore>
//http://cython.readthedocs.io/en/latest/src/userguide/wrapping_CPlusPlus.html

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double

MonoElementsIntegralCpp::MonoElementsIntegralCpp (): nodes(0),connectivity(0),lnumbering() {
  this->totalTestDofs = 0;
  this->totalUnkownDofs = 0;
  //this->numberOfVIJ = 0;
  this->hasnormal = false;
  this->totalvijcpt = 0;
  this->onlyUpper = false;

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
void MonoElementsIntegralCpp::SetNumberOfConstants(const int& n){
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
void MonoElementsIntegralCpp::SetPoints(double* pd, const int& rows, const int& columns){
    if(this->nodes) delete this->nodes;
    this->nodes = new MapMatrixDDD(pd,rows,columns);
}




void MonoElementsIntegralCpp::SetConnectivity( INT_TYPE* pd, const int& rows, const int& columns){
    if(this->connectivity) delete this->connectivity;
    this->connectivity = new MapMatrixIDD(pd,rows,columns);
}
//

void GetJackAndDet(MapMatrixDDD& valdphidxi,
                   MatrixDD3&  xcoor,
                   int& Dimensionality,
                   MatrixDDD&  Jack,
                   double&    Jdet,
                   MatrixDDD&  Jinv){

    Jack = valdphidxi*xcoor;

    if(Dimensionality == xcoor.cols()) {

      Jinv.resize(Dimensionality,xcoor.cols());

      if(Dimensionality == 3){
        //Jdet = Jack.determinant();
        MatrixDDD& m = Jack;
        Jdet = m(0, 0) * (m(1, 1) * m(2, 2) - m(2, 1) * m(1, 2)) -
               m(0, 1) * (m(1, 0) * m(2, 2) - m(1, 2) * m(2, 0)) +
               m(0, 2) * (m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0));

        const double invdet = 1 / Jdet;

        MatrixDDD& minv = Jinv; // inverse of matrix m
        minv(0, 0) = (m(1, 1) * m(2, 2) - m(2, 1) * m(1, 2)) * invdet;
        minv(0, 1) = (m(0, 2) * m(2, 1) - m(0, 1) * m(2, 2)) * invdet;
        minv(0, 2) = (m(0, 1) * m(1, 2) - m(0, 2) * m(1, 1)) * invdet;
        minv(1, 0) = (m(1, 2) * m(2, 0) - m(1, 0) * m(2, 2)) * invdet;
        minv(1, 1) = (m(0, 0) * m(2, 2) - m(0, 2) * m(2, 0)) * invdet;
        minv(1, 2) = (m(1, 0) * m(0, 2) - m(0, 0) * m(1, 2)) * invdet;
        minv(2, 0) = (m(1, 0) * m(2, 1) - m(2, 0) * m(1, 1)) * invdet;
        minv(2, 1) = (m(2, 0) * m(0, 1) - m(0, 0) * m(2, 1)) * invdet;
        minv(2, 2) = (m(0, 0) * m(1, 1) - m(1, 0) * m(0, 1)) * invdet;
        return;
      } if(Dimensionality == 2) {
        //MatrixDDD& m = Jack;
        Jdet = Jack(0, 0) * Jack(1, 1) - Jack(1, 0) * Jack(0, 1);
        const double invdet = 1 /Jdet;

        Jinv(0,0) =  Jack(1,1)* invdet;
        Jinv(0,1) = -Jack(1,0)* invdet;
        Jinv(1,0) = -Jack(0,1)* invdet;
        Jinv(1,1) =  Jack(0,0)* invdet;

      }


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

    //Jinv.compute(Jack);
    std::cerr << "Error in the calculation of the jackob" << std::endl;
    exit(0);
  };


void GetJackAndDet(MapMatrixDDD& valdphidxi,
                   MatrixDD3&  xcoor,
                   int& Dimensionality,
                   MatrixDDD&  Jack,
                   double&    Jdet,
                   QRType&  Jinv){

    Jack = valdphidxi*xcoor;


    if(Dimensionality > xcoor.cols()){
    }

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
    } else {
        std::cerr << "error in the normal " << std::endl;
        exit(1);
    }
    /*
    //# Edge in 2D
    if Jack.shape[0] == 1 and Jack.shape[1] == 2 :
            res = np.array([Jack[1,:] -Jack[0,:]],dtype =np.float)
        # surface in 3D
    elif Jack.shape[0] == 2 and Jack.shape[1] == 3 :
            res =  np.cross(Jack[0,:],Jack[1,:])
    else:
            std::cout << "Shape of Jacobian not coherent. Possible error: an elset has the same name of the considered faset" << std::endl;

        res /= np.linalg.norm(res)
        return res

        */
        };
//
void MonoElementsIntegralCpp::Integrate( WeakForm* wform, std::vector<int>& idstotreat){


    bool hasright = false;
    const int NumberOfTerms = wform->GetNumberOfTerms();
    const int NumberOfIntegrationPoints = this->iw.rows();
    const int nbcols = this->connectivity->cols();
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
    MatrixDD3 xcoor;
    xcoor.resize(nbcols,3);

    //PRINT(this->geoSpaceNumber);



    MatrixDDD ElementMatrix(maxsizelocalTestDofs,maxsizelocalUnkownDofs);

    LocalSpace& geoSpace = this->lspaces[this->geoSpaceNumber];
    const int elemdim = geoSpace.dimensionality;

    const int idstotreat_s = idstotreat.size();
    QRType Jinv;
    int n;
    int rightIndex=0;
    int leftIndex=0;
//    std::cout << "maxsizelocalTestDofs " << maxsizelocalTestDofs<< std::endl;
//    std::cout << "maxsizelocalUnkownDofs " << maxsizelocalUnkownDofs<< std::endl;
    for(int elem_counter =0; elem_counter< idstotreat_s; ++elem_counter){
        ElementMatrix =  MatrixDDD::Zero(maxsizelocalTestDofs,maxsizelocalUnkownDofs);

        n = idstotreat[elem_counter];
//        INT_TYPE fillcpt = 0;

//         the coordinates of the nodes
        for( int j = 0; j < nbcols;++j){
            xcoor.row(j) = this->nodes->row((*this->connectivity)(n,j));
        }
        // for the moment we must copy the data in a for loop, the next line will work
        // with eigen 3.3.5
        //https://stackoverflow.com/questions/40074738/eigen-extracting-submatrix-from-vector-of-indices
        //auto xcoor = this->nodes(this->connectivity.row(n),Eigen::placeholders::all);

        for(int ip = 0; ip< NumberOfIntegrationPoints;  ++ip){
//            """ we recover the jacobian matrix """
            MatrixDDD Jack ;
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
                }
                hasright = false;


                const int numberOfProds = monom.prod.size();
                for(int prodn=0; prodn< numberOfProds; ++prodn){

                    WeakTerm& term = monom.prod[prodn];

                    if(term.internalType == 0){
                        factor *= normal(term.derDegree,0);
                        continue;
                    } else if( term.internalType == 1){
                        factor *= this->constants(term.valuesIndex_,0);
                        continue;
                    } else if( term.internalType == 2){
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
                    } else if( term.internalType == 3){
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
                    } else if( term.internalType == 4){
                        LocalSpace& cs = this->lspaces[term.spaceIndex_];
                        if(term.derDegree == 1){
                            center = cs.GetBxByBz().row(term.derCoordIndex_);
                        } else {
                            center = cs.GetNxNyNz();
                        }
                        centerNumbering = this->lnumbering[term.numberingIndex_]->row(n);
                        const int ss = centerNumbering.rows();
                        vals.resize(ss,1);
                        for(int c=0; c < ss; ++c){
                            vals(c,0) = (*(this->values[term.valuesIndex_]))(centerNumbering(c),0);
                        }

                        factor *= center.dot(vals);
                        continue;
                    } else {
                        std::cout << "Cant treat term " << term.fieldName << std::endl;
                        exit(1);
                    }
                }
                if(factor == 0) continue;

                factor *= iw(ip,0);

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
                    const int loff =   this->localTestDofsOffset[leftIndex];
                    const int roff =   this->localUnkownDofsOffset[rightIndex];
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
                            const int ii = leftNumbering[i];
                            const int jj = rightNumbering[j];

                            if(this->onlyUpper && ii > jj ){
                                continue;
                            }

                            vK[totalvijcpt] = val;
                            iK[totalvijcpt] = ii;
                            jK[totalvijcpt] = jj;
                            ++totalvijcpt;

                        }
//                        std::cout << "mII[" << leftNumbering[i] <<
//                                        "," << rightNumbering[j] <<
//                                        "] = "  << ElementMatrix(localTestDofsOffset[vi]+i,localUnkownDofsOffset[vj]+j)<< std::endl;

                    }
                }
            }
        }

//        if( fillcpt){
//            SpMat A(this->totalTestDofs,this->totalUnkownDofs);
//            container cont;
//            cont.ev = &ev[0];
//            cont.ei = &ei[0];
//            cont.ej = &ej[0];
//            cont.size = fillcpt;
//
//            A.setFromTriplets(cont.begin(), cont.end());
//            A.makeCompressed ();
////            std::cout << "A" << std::endl;
//            std::cout << A ;
//            for (int k=0; k<A.outerSize(); ++k){
//                for (SpMat::InnerIterator it(A,k); it; ++it){
//                    std::cout << "m[" << it.row() << "," <<it.col() << "] = "  << it.value() << std::endl;
//                    vK[totalvijcpt] = it.value();
//                    iK[totalvijcpt] = it.row();
//                    jK[totalvijcpt] = it.col();
//                    ++totalvijcpt;
//                }
//            }
//        exit(0);
//        }

    }

};