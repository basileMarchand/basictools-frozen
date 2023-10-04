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

namespace BasicTools
{

class ElementsContainer{
    std::string elementType;
    CBasicIndexType globaloffset;
    std::shared_ptr<MapMatrixIDD > connectivity;
    std::shared_ptr<MapMatrixID1 > originalIds;
public:
    Tags tags;
    ElementsContainer(const std::string& elemtype = "" );
    CBasicIndexType GetNumberOfElements() const ;

    MAPSETGET_MatrixIDD(Connectivity,connectivity)
    MAPSETGET_MatrixID1(OriginalIds,originalIds)

    std::string GetElementType() const ;
    template<typename T>
    void AddTag(std::string& name, T& arg1);
};

struct sortbyName {
    bool operator()(std::string a, std::string b) const {
        return a < b;
    }
};

class AllElements{
public:
    std::map<std::string, ElementsContainer, sortbyName > storage;
    ElementsContainer& GetElementsOfType(const std::string& elemtype);
    //
    CBasicIndexType GetNumberOfElements() const ;
    //
    std::string ToStr();

    ElementsContainer& operator[](const std::string& key) ;

    std::map<std::string, ElementsContainer, sortbyName >::const_iterator begin() const {return this->storage.begin();}
    std::map<std::string, ElementsContainer, sortbyName >::const_iterator end() const {return this->storage.end();}
};

class UnstructuredMesh{
    std::shared_ptr<MapMatrixDDD > nodes;
    std::shared_ptr<MapMatrixID1 > originalIDNodes;
    Tags nodesTags;
    std::map<std::string,float> props;
public:
    AllElements elements;
    UnstructuredMesh();

    MAPSETGET_MatrixID1(OriginalIds,originalIDNodes)
    MAPSETGET_MatrixDDD(Nodes,nodes)

    CBasicIndexType GetNumberOfNodes() const ;
    template<typename T>
    void AddNodalTag(std::string& name, T& arg1);

    template<typename T, typename T2>
    void AddElements(std::string& elementType, T& arg1, T2& arg2);

    template<typename T>
    void AddElementTag(std::string& elementType, std::string& tagname, T& arg1);

    void Print();

    std::string ToStr();

};

} // namespace BasicTools

#include <Containers/UnstructuredMesh.hpp>

