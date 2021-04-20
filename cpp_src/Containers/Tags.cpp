//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//
#include <Containers/Tags.h>
#include <Helpers/ToString.h>

namespace BasicTools
{
    
Tag::Tag() {
    this->ids.reset(new std::shared_ptr<MapMatrixID1 >::element_type(nullptr, 0, 1));
};

void Tag::SetName(const std::string& name) {
    this->name = name;
}
std::string Tag::GetName() const {
    return this->name;
};


INT_TYPE Tag::GetSize()const {
    if (this->ids== 0) return 0;
    return this->ids->rows();
};
//**********************************
std::string Tags::ToStr() const {
    std::string res= "   Tags : ";
    for (auto const& x : this->storage) {
        res += "(" + x.first + ":" + ToString(x.second.GetSize()) +") ";
    }
    res += "\n";
    return res;
};

void Tags::AddTag(Tag& tag) {
    this->storage[tag.GetName()] = tag;
};
bool Tags::Contains(const std::string& key) const {
    //return this->storage.contains(key) ;
    return (this->storage.count(key) !=0) ;
};
Tag& Tags::operator[](const std::string& key) {
    return this->storage[key];
}

std::map<std::string,Tag>::iterator Tags::begin() {
    return this->storage.begin();
}
std::map<std::string,Tag>::iterator Tags::end() {
    return this->storage.end();
}

} // namespace BasicTools
