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
const DofAttachment& ElementSpace::GetDofAttachment(const int& dofNumber) const {
    return this->storage[dofNumber];
};

void ElementSpace::AppendDofAttachment(const char& entity, const int& entityNumber, const int& extraKey){
    storage.push_back( DofAttachment(entity,entityNumber, extraKey) ) ;
}

const MatrixDDD ElementSpace::GetValOfShapeFunctionsAt(const MatrixDDD& phiXiEta  ) const {
    return this->GetValOfShapeFunctionsAt(phiXiEta.coeff(0,0), phiXiEta.coeff(1,0),phiXiEta.coeff(2,0));
};

const MatrixDDD ElementSpace::GetValOfShapeFunctionsAt(const double& phi, const double&  xi, const double&  eta  ) const {
    return this->SFV(phi, xi, eta);
};

const MatrixDDD ElementSpace::GetValOfShapeFunctionsDerAt(const MatrixDDD& phiXiEta  ) const {
    return this->GetValOfShapeFunctionsDerAt(phiXiEta.coeff(0,0), phiXiEta.coeff(1,0),phiXiEta.coeff(2,0));
};

const MatrixDDD ElementSpace::GetValOfShapeFunctionsDerAt(const double& phi, const double&  xi, const double&  eta  ) const {
    return this->SFDV(phi, xi, eta);
};

// ***************  Space ******************
CBasicIndexType Space::GetNumberOfShapeFunctionsFor(const std::string& elementType){
    return this->storage[elementType].GetNumberOfShapeFunctions();
};
//
void Space::AddDofTo(const std::string& elementType, const char& entity, const int& entityNumber, const int& extraKey){
    this->storage[elementType].AppendDofAttachment(entity, entityNumber, extraKey);
}
//
const ElementSpace& Space::GetSpaceFor(const std::string& elementType)  {
    return this->storage[elementType];
}
//
void Space::Print(){
}
//
std::string ToStr(){
    std::string res = "Space \n";
    return res;
}
// ******************* space *************************
std::map<std::string, SpaceAtIP> EvaluateSpaceAt(const Space& sp, const SpaceIntegrationRule&  sir){
    std::map<std::string, SpaceAtIP> res;
    for(const  auto& item: sp.storage ){
        std::string key = item.first;
        res[key] =  EvaluateSpaceAt( item.second, sir.GetIR(key));
    }
    return res;
};

SpaceAtIP EvaluateSpaceAt(const ElementSpace& es , const IntegrationRule& ir ){
    SpaceAtIP res;
    const int nbp = ir.GetNumberOfPoints();
    for(int i=0; i<nbp; ++i){
        res.SFV = es.GetValOfShapeFunctionsAt(ir.p.row(i));
        res.SFDV = es.GetValOfShapeFunctionsDerAt(ir.p.row(i));
    }
    return res;
};

} // namespace BasicTools



