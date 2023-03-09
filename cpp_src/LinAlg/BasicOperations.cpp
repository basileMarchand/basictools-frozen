//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//
#include <LinAlg/BasicOperations.h>

#include <iostream>
#include <string>

namespace BasicTools {

//
void GetInv_Jacobian_Det(MapMatrixDDD& valdphidxi,
                         const Eigen::Ref<const MatrixDDD>& xcoor,
                         int& Dimensionality, MatrixDDD& Jack, double& Jdet,
                         Eigen::ColPivHouseholderQR<MatrixDDD>& Jinv) {
  Jack = valdphidxi * xcoor;

  if (Dimensionality == xcoor.cols()) {
    Jdet = Jack.determinant();

  } else {
    switch (Dimensionality) {
      case (0): {
        Jdet = 1;
        break;
      }
      case (1): {
        // we have and edge in 2D or 3D
        Jdet = Jack.norm();
        break;
      }
      case (2): {
        // TODO check
        Jdet = Jack.row(0).head<3>().cross(Jack.row(1).head<3>()).norm();
        break;
      }
      default: {
        std::cerr << "Error in the calculation of " << std::endl;
        assert(0);
      }
    }
  }

  Jinv.compute(Jack);
};

void GetNormal(const int& SpaceDim, const int& elementDim, MatrixDDD& Jack,
               MatrixD31& Normal) {
  if (elementDim == 2 && SpaceDim == 3) {
    ElementHelper::SpaceDimElementDim<3, 2>::Normal(Jack, Normal);
  } else if (elementDim == 1 && SpaceDim == 2) {
    ElementHelper::SpaceDimElementDim<2, 1>::Normal(Jack, Normal);
  } else if (elementDim == 1 && SpaceDim == 3) {
    ElementHelper::SpaceDimElementDim<3, 1>::Normal(Jack, Normal);
  } else {
    std::cerr << "error in the normal " << std::endl;
    exit(1);
  }
};
}  // namespace BasicTools
