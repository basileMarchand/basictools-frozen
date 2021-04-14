//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//
#pragma once


#include <Containers/Tags.h>
#include <LinAlg/EigenTypes.h>
#include <Helpers/SetGetMacros.h>
#include <Helpers/ToString.h>

#include <string>
#include <ostream>
#include <vector>
#include <memory>
#include <iostream>


class ElementsContainer{
    std::string elementType;
    INT_TYPE globaloffset;
    std::shared_ptr<MapMatrixIDD > connectivity;
    std::shared_ptr<MapMatrixID1 > originalIds;
public:
    Tags tags;
    ElementsContainer(const std::string& elemtype = "" ){
        this->elementType = elemtype;
        this->globaloffset = 0;
        this->connectivity.reset(new MapMatrixIDD(nullptr, 0, 1));
        this->originalIds.reset(new MapMatrixID1(nullptr, 0, 1));
    }

    INT_TYPE GetNumberOfElements() const {
        return this->connectivity->rows();
    }

    MAPSETGET_MatrixIDD(Connectivity,connectivity)
    MAPSETGET_MatrixID1(Ids,originalIds)
    std::string GetElementType() const {
        return this->elementType;
    }

    template<typename T>
    void AddTag(std::string& name, T& arg1){
        Tag ntag;
        ntag.SetName(name);
        ntag.SetIds(arg1);
        this->tags.AddTag(ntag);
    }
};

class AllElements{
public:
    std::map<std::string, ElementsContainer > storage;
    ElementsContainer& GetElementsOfType(const std::string& elemtype){
        if(this->storage.count(elemtype) == 0 ){
              this->storage[elemtype] = ElementsContainer(elemtype);
        }
        return this->storage[elemtype];
    };
    //
    INT_TYPE GetNumberOfElements(){
      INT_TYPE res = 0;
      for (auto const& x : this->storage){
         res += x.second.GetNumberOfElements();
      }
      return res;
    }
    //
    std::string ToStr(){
        std::string res;
        for (auto const& x : this->storage){
            res += "   ElementsContainer,    Type : ("+ ToString(x.first) + ","+ToString(x.second.GetNumberOfElements())+ "),";
            res += x.second.tags.ToStr() ;
        }
        return res;
    }
    //
    ElementsContainer& operator[](const std::string& key) {
        return this->storage[key];
    }

    //const ElementsContainer& operator[](const std::string& key) const {
    //    return this->storage[key];
    //}
};

class NativeUnstructuredMesh{
    std::shared_ptr<MapMatrixDDD > nodes;
    std::shared_ptr<MapMatrixID1 > originalIDNodes;
    Tags nodesTags;
    std::map<std::string,float> props;
public:
    AllElements elements;
    NativeUnstructuredMesh():
     nodes(std::make_shared<MapMatrixDDD>(nullptr, 0, 3) ),
     originalIDNodes(std::make_shared<MapMatrixID1>(nullptr,0,1) ){
    };

    MAPSETGET_MatrixID1(OriginalIds,originalIDNodes)
    MAPSETGET_MatrixDDD(Nodes,nodes)

    INT_TYPE GetNumberOfNodes() const {
        return this->nodes->rows();
    };

    template<typename T>
    void AddNodalTag(std::string& name, T& arg1){
        Tag ntag;
        ntag.SetName(name);
        ntag.SetIds(arg1);
        this->nodesTags.AddTag(ntag);
    };

    template<typename T, typename T2>
    void AddElemens(std::string& elementType, T& arg1, T2& arg2){
        ElementsContainer& ec = this->elements.GetElementsOfType(elementType);
        ec.SetConnectivity(arg1);
        ec.SetIds(arg2);
    };

    template<typename T>
    void AddElementTag(std::string& elementType, std::string& tagname, T& arg1){
        ElementsContainer& ec = this->elements.GetElementsOfType(elementType);
        ec.AddTag(tagname, arg1);
    };

    void Print(){
        std::cout << this->ToStr() ;
        std::cout.flush();
    };

    std::string ToStr(){
        std::string res= "UnstructuredMesh (cpp)\n" ;
        res += "  Number Of Nodes    : " + ToString(this->nodes->rows()) + "\n";
        res += this->nodesTags.ToStr();
        res += "  Number Of Elements : " + ToString(this->elements.GetNumberOfElements()) + "\n";
        res += this->elements.ToStr();
        return res;
    };

};

