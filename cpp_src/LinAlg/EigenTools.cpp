//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//

#include <cassert>
#include <algorithm>

#include <LinAlg/EigenTools.h>

namespace BasicTools
{
    
MatrixID1 NonZero(const MatrixBD1& input){
    const INT_TYPE& rows = input.array().count();
    MatrixID1 res(rows,1);
    INT_TYPE cpt=0;
    for(INT_TYPE i=0; i < input.cols();++i  ){
        if (input(i,0)){
            res(cpt++,0) = i;
        }
    };
    return res;
};


MatrixID1 Intersect1D(const MatrixID1& A, const MatrixID1& B){
    MatrixID1 res(std::min(A.rows(),B.rows()),1);
    INT_TYPE cptA= 0;
    INT_TYPE cptB= 0;
    INT_TYPE cpt= 0;
    while(cptA < A.rows() or cptB < B.rows() ){
       if( A(cptA,0) < B(cptB,0)){
           ++cptA;
       } else if( A(cptA,0) > B(cptB,0)){
           ++cptB;
       } else{
           res(cpt,0) = B(cptB,0);
           ++cpt;
           ++cptA;
           ++cptB;
       }
    }
    res.conservativeResize(cpt,1);
    return res;
};
//
MatrixID1 Union1D(const MatrixID1& A, const MatrixID1& B){
    MatrixID1 res(A.rows()+B.rows(),1);
    INT_TYPE cptA= 0;
    INT_TYPE cptB= 0;
    INT_TYPE cpt= 0;
    while(cptA < A.rows() or cptB < B.rows() ){
       if( A(cptA,0) < B(cptB,0)){
           res(cpt,0) = A(cptA,0);
           ++cpt;
           ++cptA;
       } else if ( A(cptA,0) > B(cptB,0)){
           res(cpt,0) = B(cptB,0);
           ++cpt;
           ++cptB;
       } else{
           res(cpt,0) = B(cptB,0);
           ++cpt;
           ++cptA;
           ++cptB;
       }
    }
    for(;cptA < A.rows(); ++ cptA ){
       res(cpt,0) = A(cptA,0);
       ++cpt;
    }
    for(;cptB < B.rows(); ++ cptB ){
       res(cpt,0) = B(cptB,0);
       ++cpt;
    }
    res.conservativeResize(cpt,1);
    return res;
}
//  from http://eigen.tuxfamily.org/dox/TopicCustomizing_NullaryExpr.html

} // namespace BasicTools
