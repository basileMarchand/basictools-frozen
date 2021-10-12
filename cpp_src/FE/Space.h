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
    CBasicIndexType GetNumberOfShapeFunctions();
    const DofAttachment& GetDofAttachment(const int& dofnumber) const ;
    void AppendDofAttachement(const char& entity, const int& entityNumber, const int& extraKey);
    const MatrixDDD GetValOfShapeFunctionsAt(const double& phi, const double&  xi, const double&  eta  ) const ;
    const MatrixDDD GetValOfShapeFunctionsAt(const MatrixDDD& phixieta) const ;
    const MatrixDDD GetValOfShapeFunctionsDerAt(const double& phi, const double&  xi, const double&  eta) const  ;
    const MatrixDDD GetValOfShapeFunctionsDerAt(const MatrixDDD& phixieta) const ;
};

class Space{
public:
    std::map<std::string, ElementSpace > storage;
public:
    CBasicIndexType GetNumberOfShapeFunctionsFor(const std::string& elemtype);
    void AddDofTo(const std::string& elemtype, const char& entity, const int& entityNumber, const int& extraKey);
    const ElementSpace& GetSpaceFor(const std::string& elemtype);
    void Print();
    std::string ToStr();
};

const Space& GetSpaceFor(const std::string& spacename);

class SpaceAtIP{
public:
    MatrixDDD SFV;
    MatrixDDD SFDV;
};

std::map<std::string, SpaceAtIP> EvaluateSpaceAt(const Space&, const SpaceIntegrationRule& );

SpaceAtIP EvaluateSpaceAt(const ElementSpace&, const IntegrationRule& );

} // namespace BasicTools
