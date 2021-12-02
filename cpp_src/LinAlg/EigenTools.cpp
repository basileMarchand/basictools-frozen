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
    const CBasicIndexType rows = static_cast<CBasicIndexType>(input.array().count());
    MatrixID1 res(rows,1);
    CBasicIndexType cpt=0;
    for(CBasicIndexType i=0; i < input.cols();++i  ){
        if (input(i,0)){
            res(cpt++,0) = i;
        }
    };
    return res;
};


MatrixID1 Intersect1D(const MatrixID1& A, const MatrixID1& B){
    MatrixID1 res(std::min(A.rows(),B.rows()),1);
    CBasicIndexType cptA = 0;
    CBasicIndexType cptB = 0;
    CBasicIndexType cpt = 0;
    while((cptA < A.rows() ) || (cptB < B.rows() ) ){
       if( A(cptA,0) < B(cptB,0)){
           ++cptA;
       } else if( A(cptA,0) > B(cptB,0)){
           ++cptB;
       } else {
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
    CBasicIndexType cptA = 0;
    CBasicIndexType cptB = 0;
    CBasicIndexType cpt = 0;
    while((cptA < A.rows() ) || ( cptB < B.rows() ) ){
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

} // namespace BasicTools
