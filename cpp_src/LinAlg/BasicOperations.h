//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//

#include <LinAlg/EigenTypes.h>

namespace BasicTools {

template<typename TM, typename TIM, typename S>
void inv22(TM& m, TIM& minv, S& det) {
  minv(0, 0) = m(1, 1);
  minv(1, 0) = -m(1, 0);
  minv(0, 1) = -m(0, 1);
  minv(1, 1) = m(0, 0);

  det = m(0, 0) * m(1, 1) - m(1, 0) * m(0, 1);

  minv *= 1. / det;
};

template<typename TM, typename TIM, typename S>
void inv33(TM& m, TIM& minv, S& det) {
  minv(0, 0) = (m(1, 1) * m(2, 2) - m(2, 1) * m(1, 2));
  minv(0, 1) = (m(0, 2) * m(2, 1) - m(0, 1) * m(2, 2));
  minv(0, 2) = (m(0, 1) * m(1, 2) - m(0, 2) * m(1, 1));
  minv(1, 0) = (m(1, 2) * m(2, 0) - m(1, 0) * m(2, 2));
  minv(1, 1) = (m(0, 0) * m(2, 2) - m(0, 2) * m(2, 0));
  minv(1, 2) = (m(1, 0) * m(0, 2) - m(0, 0) * m(1, 2));
  minv(2, 0) = (m(1, 0) * m(2, 1) - m(2, 0) * m(1, 1));
  minv(2, 1) = (m(2, 0) * m(0, 1) - m(0, 0) * m(2, 1));
  minv(2, 2) = (m(0, 0) * m(1, 1) - m(1, 0) * m(0, 1));

  det = m(0, 0) * minv(0, 0) + m(0, 1) * minv(1, 0) + m(0, 2) * minv(2, 0);

  minv *= (1. / det);
};

void GetInv_Jacobian_Det(MapMatrixDDD& valdphidxi, const Eigen::Ref<const MatrixDD3>& xcoor, int& Dimensionality, MatrixDDD& Jack, double& Jdet, Eigen::ColPivHouseholderQR<MatrixDDD>& Jinv);

void GetNormal(const int& SpaceDim, const int& elementDim, MatrixDDD& Jack, MatrixD31& Normal);

namespace ElementHelper {

template <unsigned int SpaceDim, unsigned int ElementDim>
class SpaceDimElementDim ;

template <unsigned int SpaceDim>
class SpaceDimElementDim<SpaceDim, 0> {
 public:
  template<typename JT, typename IJT>
  static void InverseJacobianDeterminant( MapMatrixDDD& valdphidxi, const Eigen::Ref<const MatrixDDD>& xcoor, JT& Jack, double& Jdet, IJT& Jinv) {
    Jdet = 1.;
    Jack(0, 0) = 1.;
    Jinv(0, 0) = 1.;
    return;
  };
  template<typename JT>
  static void Normal(JT& Jack, MatrixD31& Normal) {}
  static bool CanProvideNormal() { return false; }
  constexpr static unsigned int size = 1;
  constexpr static unsigned int spaceSize = 1;
  template<typename JT, typename IJT>
  static void ResizeContainers( JT& Jack, IJT& Jinv){
    Jack.resize(1, 1);
    Jinv.resize(1, 1);
  }
};
//-------------- 1D elements ---------------
template <>
class SpaceDimElementDim<1, 1> {
 public:
  template<typename JT, typename IJT>
  static void InverseJacobianDeterminant( MapMatrixDDD& valdphidxi, const Eigen::Ref<const MatrixDDD>& xcoor, JT& Jack, double& Jdet, IJT& Jinv) {
    Jack = valdphidxi * xcoor;
    Jdet = Jack(0, 0);
    const double invdet = 1 / Jdet;
    Jinv(0, 0) = invdet;
    return;
  };
  template<typename JT>
  static void Normal(JT& Jack, MatrixD31& Normal) {}
  static bool CanProvideNormal() { return false; }
  constexpr static unsigned int size = 1U;
  constexpr static unsigned int spaceSize = 1;
  template<typename JT, typename IJT>
  static void ResizeContainers( JT& Jack, IJT& Jinv){
    Jack.resize(1, 1);
    Jinv.resize(1, 1);
  }
};

template <>
class SpaceDimElementDim<2, 1> {
 public:
  template<typename JT, typename IJT>
  static void InverseJacobianDeterminant( MapMatrixDDD& valdphidxi, const Eigen::Ref<const MatrixDDD>& xcoor, JT& Jack, double& Jdet, IJT& Jinv) {
    // in a 2D space
    Jack = valdphidxi * xcoor;

    const double der0 = Jack(0, 0);  // dx
    const double der1 = Jack(0, 1);  // dy
    Jdet = std::sqrt(der0 * der0 + der1 * der1);
    MatrixD21 Normal;
    Normal(0) =  der1/Jdet;
    Normal(1) = -der0/Jdet;

    MatrixDDD m;
    m.resize(2, 2), m(0, 0) = der0;
    m(1, 0) = der1;
    m.col(1) = Normal;

    Jinv.noalias() = m.inverse().transpose().col(0);
    return;
  };
  template<typename JT>
  static void Normal(JT& Jack, MatrixD31& Normal) {
    const double der0 = Jack(0, 0);  // dx
    const double der1 = Jack(0, 1);  // dy
    const double Jdet = std::sqrt(der0 * der0 + der1 * der1);
    Normal(0,0) =  der1/Jdet;
    Normal(1,0) = -der0/Jdet;
    Normal(2,0) = 0;

  }
  static bool CanProvideNormal() { return true; }
  constexpr static unsigned int size = 1U;
  constexpr static unsigned int spaceSize = 2U;

  template<typename JT, typename IJT>
  static void ResizeContainers( JT& Jack, IJT& Jinv){
    Jack.resize(1, 2);
    Jinv.resize(1, 2);
  }

};

template <>
class SpaceDimElementDim<3, 1> {
 public:
  template<typename JT, typename IJT>
  static void InverseJacobianDeterminant( MapMatrixDDD& valdphidxi, const Eigen::Ref<const MatrixDDD>& xcoor, JT& Jack, double& Jdet, IJT& Jinv) {
    Jack = valdphidxi * xcoor;

    const double der0 = Jack(0, 0);  // dx
    const double der1 = Jack(0, 1);  // dy
    const double der2 = Jack(0, 2);  // dz
    Jdet = std::sqrt(der0 * der0 + der1 * der1 + der2 * der2);
    MatrixD31 N0;
    MatrixD31 N1;
    if (der0 == 0 && der1 == 0) {
      N0 << 1., 0., 0.;
      N1 << 0., 1., 0.;
    } else {
      N0 << der1, -der0, 0.;
      N1 << -der2, 0., der0;
    }
    N0 *= 1. / N0.norm();
    N1 *= 1. / N1.norm();

    MatrixD33 m;
    m.resize(3, 3), m.col(0) = Jack.row(0);
    m.col(1) = N0;
    m.col(2) = N1;

    Jinv.noalias() = m.inverse().transpose().col(0);
    return;
  };
  template<typename JT>
  static void Normal(JT& Jack, MatrixD31& Normal) {
    MatrixD31 Z1;
    Z1 << 0, 0., 1;
    Normal = Jack.row(0).template head<3>().cross(Z1);
    Normal /= Normal.norm();
  }
  static bool CanProvideNormal() { return true; }
  constexpr static unsigned int size = 1U;
  constexpr static unsigned int spaceSize = 3U;
  template<typename JT, typename IJT>
  static void ResizeContainers( JT& Jack, IJT& Jinv){
    Jack.resize(1, 3);
    Jinv.resize(1, 3);
  }
};

//-------------- 2D elements ---------------
template <>
class SpaceDimElementDim<2, 2> {
 public:
  template<typename JT, typename IJT>
  static void InverseJacobianDeterminant( MapMatrixDDD& valdphidxi, const Eigen::Ref<const MatrixDDD>& xcoor, JT& Jack, double& Jdet, IJT& Jinv) {
    Jack = valdphidxi * xcoor;
    inv22(Jack, Jinv, Jdet);
    return;
  };
  template<typename JT>
  static void Normal(JT& Jack, MatrixD31& Normal) {}
  static bool CanProvideNormal() { return false; }
  constexpr static unsigned int size = 2U;
  constexpr static unsigned int spaceSize = 2U;

  template<typename JT, typename IJT>
  static void ResizeContainers( JT& Jack, IJT& Jinv){
    Jack.resize(2, 2);
    Jinv.resize(2, 2);
  }
};

template <>
class SpaceDimElementDim<3, 2> {
 public:
  template<typename JT, typename IJT>
  static void InverseJacobianDeterminant( MapMatrixDDD& valdphidxi, const Eigen::Ref<const MatrixDDD>& xcoor, JT& Jack, double& Jdet, IJT& Jinv) {
    Jack = valdphidxi * xcoor;
    MatrixD31 Normal = Jack.row(0).template head<3>().cross(Jack.row(1).template head<3>());
    Jdet = Normal.norm();
    Normal /= Jdet;
    MatrixD33 m;
    m.row(0) = Jack.row(0);
    m.row(1) = Jack.row(1);
    m.row(2) = Normal;
    Jinv.noalias() = m.inverse().block<3, 2>(0, 0);
    return;
  };

  template<typename JT>
  static void Normal(JT& Jack, MatrixD31& Normal) {
    Normal = Jack.row(0).template head<3>().cross(Jack.row(1).template head<3>());
    Normal /= Normal.norm();
  }
  static bool CanProvideNormal() { return true; }
  constexpr static unsigned int size = 2U;
  constexpr static unsigned int spaceSize = 3U;

  template<typename JT, typename IJT>
  static void ResizeContainers( JT& Jack, IJT& Jinv){
    Jack.resize(size, 3);
    Jinv.resize(size, 3);
  }
};
//-------------- 3D elements ---------------
template <>
class SpaceDimElementDim<3, 3> {
 public:
  template<typename JT, typename IJT>
  static void InverseJacobianDeterminant( MapMatrixDDD& valdphidxi, const Eigen::Ref<const MatrixDDD>& xcoor, JT& Jack, double& Jdet, IJT& Jinv) {
    Jack = valdphidxi * xcoor;
    inv33(Jack, Jinv, Jdet);
    return;
  };
  template<typename JT>
  static void Normal(JT& Jack, MatrixD31& Normal) {}
  static bool CanProvideNormal() { return false; }
  constexpr static unsigned int size = 3U;
  constexpr static unsigned int spaceSize = 3U;
  template<typename JT, typename IJT>
  static void ResizeContainers( JT& Jack, IJT& Jinv){
    Jack.resize(3, 3);
    Jinv.resize(3, 3);
  }
};

}  // namespace ElementHelper

}  // namespace BasicTools
