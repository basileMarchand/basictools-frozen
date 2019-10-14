//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//


#include <Eigen/Dense>
#include <Eigen/Sparse>


using Eigen::Dynamic;

typedef double FLOAT_TYPE;
typedef long int INT_TYPE;


typedef Eigen::Matrix<FLOAT_TYPE, 2, 1> MatrixD21;

typedef Eigen::Matrix<FLOAT_TYPE, 3, 3> MatrixD33;
typedef Eigen::Matrix<FLOAT_TYPE, 3, 1> MatrixD31;
typedef Eigen::Matrix<FLOAT_TYPE, Dynamic, 1> MatrixDD1;
typedef Eigen::Matrix<FLOAT_TYPE, Dynamic, 3,Eigen::RowMajor> MatrixDD3;
typedef Eigen::Matrix<FLOAT_TYPE, Dynamic, Dynamic,Eigen::RowMajor> MatrixDDD;

typedef Eigen::Matrix<INT_TYPE, Dynamic, Dynamic,Eigen::RowMajor> MatrixIDD;
typedef Eigen::Matrix<INT_TYPE, Dynamic, 1> MatrixID1;
typedef Eigen::Matrix<INT_TYPE, 1, Dynamic> MatrixI1D;


typedef Eigen::SparseMatrix<FLOAT_TYPE, Eigen::ColMajor> SpMatD;


typedef Eigen::Map<MatrixDDD> MapMatrixDDD;
typedef Eigen::Map<MatrixIDD> MapMatrixIDD;
typedef Eigen::Map<MatrixDD1> MapMatrixDD1;

//typedef Eigen::ColPivHouseholderQR<MatrixDDD> QRType;
//typedef Eigen::FullPivHouseholderQR<MatrixDDD> QRType;
typedef MatrixDDD QRType;

//typedef Eigen::Map<MatrixDDD, Eigen::Aligned64> MapMatrixDDD;
//typedef Eigen::Map<MatrixIDD, Eigen::Aligned64> MapMatrixIDD;
//typedef Eigen::Map<MatrixDD1, Eigen::Aligned64> MapMatrixDD1;
