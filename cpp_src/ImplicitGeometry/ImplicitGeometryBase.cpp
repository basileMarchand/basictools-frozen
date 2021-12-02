
//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//

#include <ImplicitGeometry/ImplicitGeometryBase.h>
#include <Containers/UnstructuredMeshTools.h>

namespace BasicTools
{
void ImplicitGeometryBase::ApplyInsideOut(MatrixDD1& res)const{
    if (this->insideOut)
        res *= -1;
}
//
ImplicitGeometryBase::ImplicitGeometryBase(): insideOut(false){};
//
MatrixDD1 ImplicitGeometryBase::GetDistanceToPoints(const MatrixDDD& pos) const{
    MatrixDD1 res = GetDistanceToPointsInternal(pos);
    this->ApplyInsideOut(res);
    return res;
};

MatrixDD1 ImplicitGeometryBase::GetDistanceToMesh(UnstructuredMesh& mesh,bool cellCenter) const{
    if(cellCenter){
        MatrixDD1 res;
        res.resize(mesh.elements.GetNumberOfElements());
        CBasicIndexType cpt = 0;
        for( auto& x : mesh.elements){
            MatrixDDD centers = GetElementsCenters(mesh.GetNodesMatrix(),x.second);
            res.block(cpt,0,x.second.GetNumberOfElements(),1) = this->GetDistanceToPoints(centers);
            cpt += x.second.GetNumberOfElements();
        }
        return res;
    } else {
        return this->GetDistanceToPoints(mesh.GetNodesMatrix());
    }
};

void ImplicitGeometryBase::SetInsideOutOn(){
    this->insideOut = true;
}
void ImplicitGeometryBase::SetInsideOutOff(){
    this->insideOut = false;
}
void ImplicitGeometryBase::SetInsideOut(const bool& val){
    this->insideOut = val;
}
bool ImplicitGeometryBase::GetInsideOut()const{
    return this->insideOut;
}

} // namespace BasicTools
