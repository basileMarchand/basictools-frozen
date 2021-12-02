//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//

#include <FE/Space.h>

namespace BasicTools
{

// ************ ElementSpace ****************
CBasicIndexType ElementSpace::GetNumberOfShapeFunctions(){
    return static_cast<CBasicIndexType>(this->storage.size());
}
//
const DofAttachment& ElementSpace::GetDofAttachment(const int& dofnumber) const {
    return this->storage[dofnumber];
};

void ElementSpace::AppendDofAttachement(const char& entity, const int& entityNumber, const int& extraKey){
    storage.push_back( DofAttachment(entity,entityNumber, extraKey) ) ;
}
// ***************  Space ******************
CBasicIndexType Space::GetNumberOfShapeFunctionsFor(const std::string& elemtype){
    return this->storage[elemtype].GetNumberOfShapeFunctions();
};
//
void Space::AddDofTo(const std::string& elemtype, const char& entity, const int& entityNumber, const int& extraKey){
    this->storage[elemtype].AppendDofAttachement(entity, entityNumber, extraKey);
}
//
ElementSpace& Space::GetSpaceFor(const std::string& elemtype){
    return this->storage[elemtype];
}
//
void Space::Print(){
}
//
std::string ToStr(){
    std::string res = "Space \n";
    return res;
}

} // namespace BasicTools