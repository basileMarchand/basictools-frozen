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
    CBasicIndexType GetNumberOfShapeFunctions();
    const DofAttachment& GetDofAttachment(const int& dofnumber) const ;
    void AppendDofAttachement(const char& entity, const int& entityNumber, const int& extraKey);
};

class Space{
    std::map<std::string, ElementSpace > storage;
public:
    CBasicIndexType GetNumberOfShapeFunctionsFor(const std::string& elemtype);
    void AddDofTo(const std::string& elemtype, const char& entity, const int& entityNumber, const int& extraKey);
    ElementSpace& GetSpaceFor(const std::string& elemtype);
    void Print();
    std::string ToStr();
};

} // namespace BasicTools
