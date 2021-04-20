//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//
#pragma once

#include <Helpers/SetGetMacros.h>
#include <LinAlg/EigenTypes.h>
#include <string>
#include <memory>

namespace BasicTools
{
    
class Tag {
    std::string name;
    std::shared_ptr<MapMatrixID1 >  ids;
public:
    Tag();
    MAPSETGET_MatrixID1(Ids,ids)
    void SetName(const std::string& name);
    std::string GetName() const ;
    INT_TYPE GetSize()const ;

    template<typename T>
    void SetIds(T& ids);
};

class Tags {
    std::map<std::string,Tag> storage;
public:
    void AddTag(Tag& tag);
    bool Contains(const std::string& key) const ;
    Tag& operator[](const std::string& key);
    std::string ToStr() const ;
    std::map<std::string,Tag>::iterator begin();
    std::map<std::string,Tag>::iterator end();
};

} // namespace BasicTools

#include <Containers/Tags.hpp>
