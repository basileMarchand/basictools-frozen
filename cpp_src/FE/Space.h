//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//
#pragma once

#include <string>
#include <ostream>
#include <vector>
#include <memory>

#include <LinAlg/EigenTypes.h>
#include <FE/IntegrationRule.h>

namespace BasicTools
{

struct DofAttachment {
    DofAttachment(const char& entity, const int& entityNumber, const int& extraKey): entity(entity), entityNumber(entityNumber), extraKey(extraKey){};
    char entity;
    int  entityNumber;
    int  extraKey;
};

struct ElementSpace{
    std::vector<DofAttachment > storage;
    ElementSpace(){};

    MatrixDDD (*SFV)(const double&,const double&,const double& );
    MatrixDDD (*SFDV)(const double&,const double&,const double& );
    CBasicIndexType GetNumberOfShapeFunctions() const;
    const DofAttachment& GetDofAttachment(const int& dofNumber) const ;
    void AppendDofAttachment(const char& entity, const int& entityNumber, const int& extraKey);
    const MatrixDDD GetValOfShapeFunctionsAt(const double& phi, const double&  xi, const double&  eta  ) const ;
    const MatrixDDD GetValOfShapeFunctionsAt(const MatrixDDD& phiXiEta) const ;
    const MatrixDDD GetValOfShapeFunctionsDerAt(const double& phi, const double&  xi, const double&  eta) const  ;
    const MatrixDDD GetValOfShapeFunctionsDerAt(const MatrixDDD& phiXiEta) const ;
};

class Space{
public:
    std::map<std::string, ElementSpace > storage;
public:
    CBasicIndexType GetNumberOfShapeFunctionsFor(const std::string& elementType);
    void AddDofTo(const std::string& elementType, const char& entity, const int& entityNumber, const int& extraKey);
    const ElementSpace& GetSpaceFor(const std::string& elementType) const ;
    void Print();
    std::string ToStr();
};

class SpaceAtIP{
public:
    std::vector<MatrixDDD> SFV;
    std::vector<MatrixDDD> SFDV;
};

std::map<std::string, SpaceAtIP> EvaluateSpaceAt(const Space&, const SpaceIntegrationRule& );

SpaceAtIP EvaluateSpaceAt(const ElementSpace&, const IntegrationRule& );

} // namespace BasicTools
