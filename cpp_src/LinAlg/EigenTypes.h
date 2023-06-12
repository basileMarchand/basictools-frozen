//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//

#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace BasicTools
{

#ifdef __linux__
    using CBasicIndexType = Eigen::Index;
#else
    using CBasicIndexType = long;
#endif

using CBasicFloatType = double;

typedef Eigen::Matrix<CBasicFloatType, 2, 1> MatrixD21;

typedef Eigen::Matrix<CBasicFloatType, 3, 3> MatrixD33;
typedef Eigen::Matrix<CBasicFloatType, 3, 1> MatrixD31;
typedef Eigen::Matrix<CBasicFloatType, 3, 1> MatrixD13;
typedef Eigen::Matrix<CBasicFloatType, Eigen::Dynamic, 1> MatrixDD1;
typedef Eigen::Matrix<CBasicFloatType, 1, Eigen::Dynamic> MatrixD1D;

typedef Eigen::Matrix<CBasicFloatType, Eigen::Dynamic, 3,Eigen::RowMajor> MatrixDD3;
typedef Eigen::Matrix<CBasicFloatType, Eigen::Dynamic, Eigen::Dynamic,Eigen::RowMajor> MatrixDDD;

typedef Eigen::Matrix<CBasicIndexType, Eigen::Dynamic, Eigen::Dynamic,Eigen::RowMajor> MatrixIDD;
typedef Eigen::Matrix<CBasicIndexType, Eigen::Dynamic, 1> MatrixID1;
typedef Eigen::Matrix<CBasicIndexType, 1, Eigen::Dynamic> MatrixI1D;

typedef Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic,Eigen::RowMajor> MatrixBDD;
typedef Eigen::Matrix<bool, Eigen::Dynamic, 1> MatrixBD1;
typedef Eigen::Matrix<bool, 1, Eigen::Dynamic> MatrixB1D;

// Arrays ****************************

typedef Eigen::Array<CBasicFloatType, Eigen::Dynamic, Eigen::Dynamic> ArrayDDD;
typedef Eigen::Array<CBasicFloatType, Eigen::Dynamic, 1>              ArrayDD1;

typedef Eigen::Array<CBasicIndexType, Eigen::Dynamic, Eigen::Dynamic>   ArrayIDD;
typedef Eigen::Array<CBasicIndexType, Eigen::Dynamic, 1>                ArrayID1;

typedef Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic>   ArrayBDD;
typedef Eigen::Array<bool, Eigen::Dynamic, 1>                ArrayBD1;

// SpMatrix ****************************

typedef Eigen::SparseMatrix<CBasicFloatType, Eigen::ColMajor> SpMatD;
typedef Eigen::SparseMatrix<CBasicFloatType, Eigen::RowMajor> SpMatDR;

// Maps ****************************

typedef Eigen::Map<MatrixDDD> MapMatrixDDD;
typedef Eigen::Map<MatrixIDD> MapMatrixIDD;

typedef Eigen::Map<MatrixDD1> MapMatrixDD1;
typedef Eigen::Map<MatrixID1> MapMatrixID1;

// Maps of const

typedef Eigen::Map<const MatrixDDD> MapConstMatrixDDD;
typedef Eigen::Map<const MatrixIDD> MapConstMatrixIDD;

typedef Eigen::Map<const MatrixDD1> MapConstMatrixDD1;
typedef Eigen::Map<const MatrixID1> MapConstMatrixID1;

typedef MatrixDDD QRType;

} // BasicTools namespace

