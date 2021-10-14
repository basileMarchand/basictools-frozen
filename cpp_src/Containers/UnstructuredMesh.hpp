//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//

namespace BasicTools
{
    
template<typename T>
void ElementsContainer::AddTag(std::string& name, T& arg1){
    Tag ntag;
    ntag.SetName(name);
    ntag.SetIds(arg1);
    this->tags.AddTag(ntag);
}

//////////// UnstructuredMesh

template<typename T>
void UnstructuredMesh::AddNodalTag(std::string& name, T& arg1){
    Tag ntag;
    ntag.SetName(name);
    ntag.SetIds(arg1);
    this->nodesTags.AddTag(ntag);
};

template<typename T, typename T2>
void UnstructuredMesh::AddElemens(std::string& elementType, T& arg1, T2& arg2){
    ElementsContainer& ec = this->elements.GetElementsOfType(elementType);
    ec.SetConnectivity(arg1);
    ec.SetOriginalIds(arg2);
};

template<typename T>
void UnstructuredMesh::AddElementTag(std::string& elementType, std::string& tagname, T& arg1){
    ElementsContainer& ec = this->elements.GetElementsOfType(elementType);
    ec.AddTag(tagname, arg1);
};

} // namespace BasicTools
