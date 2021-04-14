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

struct DofAttachment {
    DofAttachment(const char& entity, const int& entityNumber, const int& extraKey): entity(entity), entityNumber(entityNumber), extraKey(extraKey){};
    char entity;
    int  entityNumber;
    int  extraKey;


};

struct ElementSpace{
    std::vector<DofAttachment > storage;
    ElementSpace(){}
    int GetNumberOfShapeFunctions(){
        return this->storage.size();
    }
    const DofAttachment& GetDofAttachment(const int& dofnumber) const {
        return this->storage[dofnumber];
    };
    void AppendDofAttachement(const char& entity, const int& entityNumber, const int& extraKey){
        //std::cout << "AppendDofAttachement  "  << entity << " " << entityNumber << " " << extraKey << std::endl;
        storage.push_back( DofAttachment(entity,entityNumber, extraKey) ) ;
    }
};

class Space{
    std::map<std::string, ElementSpace > storage;
public:
    int GetNumberOfShapeFunctionsFor(const std::string& elemtype){
        return this->storage[elemtype].GetNumberOfShapeFunctions();
    };
    void AddDofTo(const std::string& elemtype, const char& entity, const int& entityNumber, const int& extraKey){
        this->storage[elemtype].AppendDofAttachement(entity, entityNumber, extraKey);
    }
    ElementSpace& GetSpaceFor(const std::string& elemtype){
        return this->storage[elemtype];
    }

    void Print(){
    }
    std::string ToStr(){
        std::string res = "Space \n";
        return res;
    }
};

