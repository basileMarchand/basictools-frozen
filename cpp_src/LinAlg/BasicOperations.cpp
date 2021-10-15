//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//
#include <iostream>
#include <string>

#include <LinAlg/BasicOperations.h>

namespace BasicTools
{

void inv22(MatrixDDD& m,MatrixDDD &minv,double& det ){

        // in 2D space
        det = m(0, 0) * m(1, 1) - m(1, 0) * m(0, 1);
        const double invdet = 1. /det;

        minv(0,0) =  m(1,1)* invdet;
        minv(1,0) = -m(1,0)* invdet;
        minv(0,1) = -m(0,1)* invdet;
        minv(1,1) =  m(0,0)* invdet;
}; 

void inv33(MatrixDDD& m,MatrixDDD &minv,double& det ){
        det = m(0, 0) * (m(1, 1) * m(2, 2) - m(2, 1) * m(1, 2)) -
               m(0, 1) * (m(1, 0) * m(2, 2) - m(1, 2) * m(2, 0)) +
               m(0, 2) * (m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0));

        //const double invdet = 1. / det;
        minv(0, 0) = (m(1, 1) * m(2, 2) - m(2, 1) * m(1, 2));// * invdet;
        minv(0, 1) = (m(0, 2) * m(2, 1) - m(0, 1) * m(2, 2));// * invdet;
        minv(0, 2) = (m(0, 1) * m(1, 2) - m(0, 2) * m(1, 1));//* invdet;
        minv(1, 0) = (m(1, 2) * m(2, 0) - m(1, 0) * m(2, 2));// * invdet;
        minv(1, 1) = (m(0, 0) * m(2, 2) - m(0, 2) * m(2, 0));// * invdet;
        minv(1, 2) = (m(1, 0) * m(0, 2) - m(0, 0) * m(1, 2));// * invdet;
        minv(2, 0) = (m(1, 0) * m(2, 1) - m(2, 0) * m(1, 1));// * invdet;
        minv(2, 1) = (m(2, 0) * m(0, 1) - m(0, 0) * m(2, 1));// * invdet;
        minv(2, 2) = (m(0, 0) * m(1, 1) - m(1, 0) * m(0, 1));// * invdet;
        minv *=  (1./det);
};

void Inv_Jacobian_Det(MapMatrixDDD& valdphidxi,
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
            inv22(Jack,Jinv,Jdet);
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
void Inv_Jacobian_Det(MapMatrixDDD& valdphidxi,
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
} // namespace BasicTools
