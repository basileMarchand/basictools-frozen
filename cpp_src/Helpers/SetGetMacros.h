//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//
#pragma once

#define MAPSETGETTYPEII(Name,attribute,dtype)                                   \
    template<typename T>                                                        \
    void Set##Name(Eigen::Map<T>  &arg1){                                       \
        this->attribute.reset(new dtype(arg1.data(), arg1.rows(), arg1.cols()));\
    };                                                                          \
    template<typename T>                                                        \
    void Set##Name(std::shared_ptr<T>  &arg1){                                  \
        this->attribute = arg1;                                                 \
    };                                                                          \
    std::shared_ptr<dtype>& Get ## Name ## Shared(){                            \
        return this->attribute;                                                 \
    };                                                                          \
    dtype& Get ## Name ## Matrix(){                                             \
        return *(this->attribute.get());                                        \
    };                                                                          \
    const dtype& Get ## Name ## Matrix() const {                                \
        return *(this->attribute.get());                                        \
    };

#define MAPSETGET_MatrixID1(Name,attribute) MAPSETGETTYPEII(Name,attribute,MapMatrixID1)
#define MAPSETGET_MatrixIDD(Name,attribute) MAPSETGETTYPEII(Name,attribute,MapMatrixIDD)
#define MAPSETGET_MatrixDDD(Name,attribute) MAPSETGETTYPEII(Name,attribute,MapMatrixDDD)
