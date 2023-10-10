//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//
#pragma once

#define MAPSETGETTYPEII(Name,attribute,dtype,btype)                             \
    template<typename T>                                                        \
    void Set##Name(Eigen::Map<T>  &arg1){                                       \
        this->attribute.reset(new dtype(arg1.data(), arg1.rows(), arg1.cols()));\
    };                                                                          \
    void Set##Name(std::shared_ptr<btype>  &arg1){                              \
        this->attribute.reset(new dtype(arg1->data(), arg1->rows(), arg1->cols()));\
    };                                                                          \
    template<typename T>                                                        \
    void Set##Name(std::shared_ptr<T>  &arg1){                                  \
        this->attribute = arg1;                                                 \
    };                                                                          \
    std::shared_ptr<dtype> Get ## Name ## Shared(){                             \
        return this->attribute;                                                 \
    };                                                                          \
    dtype Get ## Name ## Matrix(){                                              \
        return *(this->attribute.get());                                        \
    };                                                                          \
    const dtype Get ## Name ## Matrix() const {                                 \
        return *(this->attribute.get());                                        \
    };

#define MAPSETGET_MatrixID1(Name,attribute) MAPSETGETTYPEII(Name,attribute,MapMatrixID1,MatrixID1)
#define MAPSETGET_MatrixIDD(Name,attribute) MAPSETGETTYPEII(Name,attribute,MapMatrixIDD,MatrixIDD)
#define MAPSETGET_MatrixDDD(Name,attribute) MAPSETGETTYPEII(Name,attribute,MapMatrixDDD,MatrixDDD)


#define MACRO_SetGet_EIGEN(Name,attribute,dtype)                                \
                                                                                \
    void Set##Name(const dtype &  &arg1){                                       \
        this->attribute.reset(new dtype(arg1.data(), arg1.rows(), arg1.cols()));\
    };                                                                          \
    dtype& Get ## Name ## Matrix(){                                             \
        return *(this->attribute);                                              \
    };                                                                          \
    const dtype& Get ## Name ## Matrix() const {                                \
        return *(this->attribute);                                              \
    };

#define SETGET_MatrixID1(Name,attribute) MACRO_SetGet_EIGEN(Name,attribute,MarixID1)
#define SETGET_MatrixIDD(Name,attribute) MACRO_SetGet_EIGEN(Name,attribute,MatrixIDD)
#define SETGET_MatrixDDD(Name,attribute) MACRO_SetGet_EIGEN(Name,attribute,MatrixDDD)
