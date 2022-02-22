//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//
#pragma once
#include <memory>

#include <LinAlg/EigenTypes.h>

namespace BasicTools {


class IntegrationRule {
public:
    MatrixDDD p;
    MatrixDDD w;
    inline CBasicIndexType GetNumberOfPoints() const {
        return static_cast<CBasicIndexType>(this->p.rows());
    }
};

class SpaceIntegrationRule {
public:
    std::map<std::string,IntegrationRule> storage;
    const IntegrationRule& GetIR(const std::string& key ) const {
        return this->storage.at(key);
    }
};


};// BasicTools namespace

