//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//

#include <LinAlg/EigenTypes.h>

namespace BasicTools
{

void inv22(MatrixDDD& m,MatrixDDD &minv,double& det );

void inv33(MatrixDDD& m,MatrixDDD &minv,double& det );

void GetInv_Jacobian_Det(MapMatrixDDD& valdphidxi, MatrixDDD& xcoor, int& Dimensionality, MatrixDDD& Jack, double& Jdet, MatrixDDD& Jinv);

void GetInv_Jacobian_Det(MapMatrixDDD& valdphidxi, MatrixDD3& xcoor, int& Dimensionality, MatrixDDD& Jack, double& Jdet, Eigen::ColPivHouseholderQR<MatrixDDD>& Jinv);

void GetNormal(const int&  SpaceDim, const int& elementDim,MatrixDDD& Jack,MatrixD31& Normal);

} // namespace BasicTools
