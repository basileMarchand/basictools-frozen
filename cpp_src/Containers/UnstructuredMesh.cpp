//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//

#include <Containers/UnstructuredMesh.h>

namespace BasicTools
{
    
//////////// ElementsContainer
ElementsContainer::ElementsContainer(const std::string& elemtype ){
    this->elementType = elemtype;
    this->globaloffset = 0;
    this->connectivity.reset(new MapMatrixIDD(nullptr, 0, 1));
    this->originalIds.reset(new MapMatrixID1(nullptr, 0, 1));
}

INT_TYPE ElementsContainer::GetNumberOfElements() const {
    return this->connectivity->rows();
}

std::string ElementsContainer::GetElementType() const {
    return this->elementType;
}

//////////// AllElements

ElementsContainer& AllElements::GetElementsOfType(const std::string& elemtype){
    if(this->storage.count(elemtype) == 0 ){
            this->storage[elemtype] = ElementsContainer(elemtype);
    }
    return this->storage[elemtype];
};
//
INT_TYPE AllElements::GetNumberOfElements() const {
    INT_TYPE res = 0;
    for (auto const& x : this->storage){
        res += x.second.GetNumberOfElements();
    }
    return res;
}
//
std::string AllElements::ToStr(){
    std::string res;
    for (auto const& x : this->storage){
        res += "   ElementsContainer,    Type : ("+ ToString(x.first) + ","+ToString(x.second.GetNumberOfElements())+ "),";
        res += x.second.tags.ToStr() ;
    }
    return res;
}
//
ElementsContainer& AllElements::operator[](const std::string& key) {
    return this->storage[key];
}
//////////// NativeUnstructuredMesh

UnstructuredMesh::UnstructuredMesh():
     nodes(std::make_shared<MapMatrixDDD>(nullptr, 0, 3) ),
     originalIDNodes(std::make_shared<MapMatrixID1>(nullptr,0,1) ){
};

INT_TYPE UnstructuredMesh::GetNumberOfNodes() const {
    return this->nodes->rows();
};

void UnstructuredMesh::Print(){
    std::cout << this->ToStr() ;
    std::cout.flush();
};

std::string UnstructuredMesh::ToStr(){
    std::string res= "UnstructuredMesh (cpp)\n" ;
    res += "  Number Of Nodes    : " + ToString(this->nodes->rows()) + "\n";
    res += this->nodesTags.ToStr();
    res += "  Number Of Elements : " + ToString(this->elements.GetNumberOfElements()) + "\n";
    res += this->elements.ToStr();
    return res;
};



} // namespace BasicTools
