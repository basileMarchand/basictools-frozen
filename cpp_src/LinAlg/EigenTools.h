//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//
#pragma once

#include <cassert>

#include <LinAlg/EigenTypes.h>

#define max(a,b) ((a)>=(b)?(a):(b))
#define min(a,b) ((a)<=(b)?(a):(b))

namespace BasicTools
{
    
MatrixID1 NonZero(const MatrixBD1& input);

MatrixID1 Intersect1D(const MatrixID1& A, const MatrixID1& B);

MatrixID1 Union1D(const MatrixID1& A, const MatrixID1& B);

//
//  from http://eigen.tuxfamily.org/dox/TopicCustomizing_NullaryExpr.html
// usage:
//  Eigen::MatrixXi A = Eigen::MatrixXi::Random(4,4);
//  Array3i ri(1,2,1);
//  ArrayXi ci(6); ci << 3,2,1,0,0,2;
//  Eigen::MatrixXi B = indexingij(A, ri, ci); // numpy like B = A[ri,:][:,ci]
//  Eigen::MatrixXi B = indexingi(A, ri, 0); // numpy like B = A[ri,0]
//  Eigen::MatrixXi B = indexingj(A, 1, ci); // numpy like B = A[1,ci]
//***********************  Indexing A(vec,vec)********************
template<class ArgType, class RowIndexType, class ColIndexType>
class indexing_functorij {
  const ArgType &m_arg;
  const RowIndexType &m_rowIndices;
  const ColIndexType &m_colIndices;
public:
  typedef Eigen::Matrix<typename ArgType::Scalar,
                 RowIndexType::SizeAtCompileTime,
                 ColIndexType::SizeAtCompileTime,
                 ArgType::Flags&Eigen::RowMajorBit?Eigen::RowMajor:Eigen::ColMajor,
                 RowIndexType::MaxSizeAtCompileTime,
                 ColIndexType::MaxSizeAtCompileTime> MatrixType;

  indexing_functorij(const ArgType& arg, const RowIndexType& row_indices, const ColIndexType& col_indices)
    : m_arg(arg), m_rowIndices(row_indices), m_colIndices(col_indices)
  {}

  const typename ArgType::Scalar& operator() (Eigen::Index row, Eigen::Index col) const {
    return m_arg(m_rowIndices[row], m_colIndices[col]);
  }
};

template <class ArgType, class RowIndexType, class ColIndexType>
Eigen::CwiseNullaryOp<indexing_functorij<ArgType,RowIndexType,ColIndexType>, typename indexing_functorij<ArgType,RowIndexType,ColIndexType>::MatrixType>
indexingij(const Eigen::MatrixBase<ArgType>& arg, const RowIndexType& row_indices, const ColIndexType& col_indices)
{
  typedef indexing_functorij<ArgType,RowIndexType,ColIndexType> Func;
  typedef typename Func::MatrixType MatrixType;
  return MatrixType::NullaryExpr(row_indices.size(), col_indices.size(), Func(arg.derived(), row_indices, col_indices));
}
//***********************  Indexing A(vec,int)********************
template<class ArgType, class RowIndexType>
class indexing_functori {
  const ArgType &m_arg;
  const RowIndexType &m_rowIndices;
  const Eigen::Index m_colIndex;

public:
  typedef Eigen::Matrix<typename ArgType::Scalar,
                 RowIndexType::SizeAtCompileTime,
                 1,
                 ArgType::Flags&Eigen::RowMajorBit?Eigen::RowMajor:Eigen::ColMajor,
                 RowIndexType::MaxSizeAtCompileTime,
                 1> MatrixType;

  indexing_functori(const ArgType& arg, const RowIndexType& row_indices, const Eigen::Index& col_index)
    : m_arg(arg), m_rowIndices(row_indices), m_colIndex(col_index)
  {}

  const typename ArgType::Scalar& operator() (Eigen::Index row, Eigen::Index col) const {
    assert(col == 0);
    return m_arg(m_rowIndices(row), m_colIndex);
  }
};

template <class ArgType, class RowIndexType>
Eigen::CwiseNullaryOp<indexing_functori<ArgType,RowIndexType>, typename indexing_functori<ArgType,RowIndexType>::MatrixType>
indexingi(const Eigen::MatrixBase<ArgType>& arg, const RowIndexType& row_indices, const Eigen::Index& col_index)
{
  typedef indexing_functori<ArgType,RowIndexType> Func;
  typedef typename Func::MatrixType MatrixType;
  return MatrixType::NullaryExpr(row_indices.size(), 1 , Func(arg.derived(), row_indices, col_index));
}

//***********************  Indexing A(int,vec)********************
template<class ArgType, class ColIndexType>
class indexing_functorj {
  const ArgType &m_arg;
  const Eigen::Index m_rowIndex;
  const ColIndexType &m_colIndices;
public:
  typedef Eigen::Matrix<typename ArgType::Scalar,
                 1,
                 ColIndexType::SizeAtCompileTime,
                 ArgType::Flags&Eigen::RowMajorBit?Eigen::RowMajor:Eigen::ColMajor,
                 1,
                 ColIndexType::MaxSizeAtCompileTime> MatrixType;

  indexing_functorj(const ArgType& arg, const Eigen::Index& row_index, const ColIndexType& col_indices)
    : m_arg(arg), m_rowIndex(row_index), m_colIndices(col_indices)
  {}

  const typename ArgType::Scalar& operator() (Eigen::Index row, Eigen::Index col) const {
    return m_arg(m_rowIndex, m_colIndices[col]);
  }
};

template <class ArgType, class ColIndexType>
Eigen::CwiseNullaryOp<indexing_functorj<ArgType,ColIndexType>, typename indexing_functorj<ArgType,ColIndexType>::MatrixType>
indexingj(const Eigen::MatrixBase<ArgType>& arg, const Eigen::Index& row_index, const ColIndexType& col_indices)
{
  typedef indexing_functorj<ArgType,ColIndexType> Func;
  typedef typename Func::MatrixType MatrixType;
  return MatrixType::NullaryExpr(1, col_indices.size(), Func(arg.derived(), row_index, col_indices));
}



// ***********************  Indexing A(mat)********************
template<class ArgType, class RowColIndexType>
class indexing_functorm {
  const ArgType &m_arg;
  const RowColIndexType &m_rowcolIndices;
public:
  typedef Eigen::Matrix<typename ArgType::Scalar,
                 RowColIndexType::SizeAtCompileTime,
                 RowColIndexType::SizeAtCompileTime,
                 ArgType::Flags&Eigen::RowMajorBit?Eigen::RowMajor:Eigen::ColMajor,
                 RowColIndexType::MaxSizeAtCompileTime,
                 RowColIndexType::MaxSizeAtCompileTime> MatrixType;

  indexing_functorm(const ArgType& arg, const RowColIndexType& rowcol_indices)
    : m_arg(arg), m_rowcolIndices(rowcol_indices)
  {}

  const typename ArgType::Scalar& operator() (Eigen::Index row, Eigen::Index col) const {
    return m_arg(m_rowcolIndices(row,col), 0);
  }
};



template <class ArgType, class RowColIndexType>
Eigen::CwiseNullaryOp<indexing_functorm<ArgType,RowColIndexType>, typename indexing_functorm<ArgType,RowColIndexType>::MatrixType>
indexingm(const Eigen::MatrixBase<ArgType>& arg, const RowColIndexType& rowcol_indices)
{
  typedef indexing_functorm<ArgType,RowColIndexType> Func;
  typedef typename Func::MatrixType MatrixType;
  return MatrixType::NullaryExpr(rowcol_indices.rows(), rowcol_indices.cols(), Func(arg.derived(), rowcol_indices));
}

//***********************   ********************

} // namespace BasicTools
