//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//

#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace BasicTools
{

typedef double FLOAT_TYPE;
typedef long int INT_TYPE;

// Matrix ****************************

typedef Eigen::Matrix<FLOAT_TYPE, 2, 1> MatrixD21;

typedef Eigen::Matrix<FLOAT_TYPE, 3, 3> MatrixD33;
typedef Eigen::Matrix<FLOAT_TYPE, 3, 1> MatrixD31;
typedef Eigen::Matrix<FLOAT_TYPE, Eigen::Dynamic, 1> MatrixDD1;
typedef Eigen::Matrix<FLOAT_TYPE, Eigen::Dynamic, 3,Eigen::RowMajor> MatrixDD3;
typedef Eigen::Matrix<FLOAT_TYPE, Eigen::Dynamic, Eigen::Dynamic,Eigen::RowMajor> MatrixDDD;

typedef Eigen::Matrix<INT_TYPE, Eigen::Dynamic, Eigen::Dynamic,Eigen::RowMajor> MatrixIDD;
typedef Eigen::Matrix<INT_TYPE, Eigen::Dynamic, 1> MatrixID1;
typedef Eigen::Matrix<INT_TYPE, 1, Eigen::Dynamic> MatrixI1D;

typedef Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic,Eigen::RowMajor> MatrixBDD;
typedef Eigen::Matrix<bool, Eigen::Dynamic, 1> MatrixBD1;
typedef Eigen::Matrix<bool, 1, Eigen::Dynamic> MatrixB1D;

// Arrays ****************************

typedef Eigen::Array<FLOAT_TYPE, Eigen::Dynamic, Eigen::Dynamic> ArrayDDD;
typedef Eigen::Array<FLOAT_TYPE, Eigen::Dynamic, 1>              ArrayDD1;

typedef Eigen::Array<INT_TYPE, Eigen::Dynamic, Eigen::Dynamic>   ArrayIDD;
typedef Eigen::Array<INT_TYPE, Eigen::Dynamic, 1>                ArrayID1;

typedef Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic>   ArrayBDD;
typedef Eigen::Array<bool, Eigen::Dynamic, 1>                ArrayBD1;

// SpMatrix ****************************

typedef Eigen::SparseMatrix<FLOAT_TYPE, Eigen::ColMajor> SpMatD;
typedef Eigen::SparseMatrix<FLOAT_TYPE, Eigen::RowMajor> SpMatDR;

// Maps ****************************

typedef Eigen::Map<MatrixDDD> MapMatrixDDD;
typedef Eigen::Map<MatrixIDD> MapMatrixIDD;

typedef Eigen::Map<MatrixDD1> MapMatrixDD1;
typedef Eigen::Map<MatrixID1> MapMatrixID1;

typedef MatrixDDD QRType;

} // BasicTools namespace

