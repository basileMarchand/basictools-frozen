//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//

namespace BasicTools
{
 
template<typename T>
void Tag::SetIds(T& ids){
    this->ids.reset(new std::shared_ptr<MapMatrixID1 >::element_type(ids.data(), ids.rows(), ids.cols()));
};

} // namespace BasicTools
