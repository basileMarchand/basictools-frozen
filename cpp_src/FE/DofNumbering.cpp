//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//

#include <FE/DofNumbering.h>
#include <algorithm>

namespace BasicTools
{
const std::string Point = "P";
const std::string Global = "G";
const std::string IntegrationPoint = "I";

//****************** DofKey ***************

DofKey::DofKey(const std::string s,const CBasicIndexType& id0, const int& id1): s(s), id1(id1){
    this->id0.resize(1,1);
    this->id0(0,0)= id0;
};

bool DofKey::operator<( DofKey const& b) const {
   if(this->s < b.s){
       return true;
   } else if (this->s == b.s ) {
       if (this->id1 < b.id1 ){
           return true;
       } else if (this->id1 == b.id1 ) {

           if( this->id0.rows() < b.id0.rows() ) return true;
           //if( this->id0.rows() >= b.id0.rows() ) return false;
           if( this->id0.rows() == b.id0.rows() ) {
               const CBasicIndexType& m = static_cast<CBasicIndexType>(this->id0.rows());
               for(CBasicIndexType i =0;  i < m ; ++i){
                   if( this->id0(i,0) < b.id0(i,0) ) return true;
                   if( this->id0(i,0) ==  b.id0(i,0) ) continue;
                   return false;
               }
           }
       }
   }
   return false;
};
void DofKey::SortId0(){
    std::sort(this->id0.data(),this->id0.data()+this->id0.size());
};
void DofKey::Print() const {
    std::cout << "DofKey "  << this->s << " " << this->id0.transpose() << "  : " <<  this->id1 <<  std::endl;
}
//
bool DofKey::IsPointDof()const{
    return Point == this->s;
}
CBasicIndexType DofKey::GetPoit()const {
    return this->id0(0,0);
}
//------------------------------------------------------------
DofKey GetPointDofKey(const CBasicIndexType&pid ) {
    return DofKey(Point,pid,0);
}

//-------------------  DofNumberingCpp -----------------------------------------
DofNumbering::DofNumbering(){
    this->size=0;
    this->fromConnectivity=false;
    this->discontinuous=false;
    this->dofToPointComputed=false;
    this->dofToCellComputed=false;
};
//
CBasicIndexType DofNumbering::GetSize() const {
    return this->size;
};
//
bool DofNumbering::GetFromConnectivity() const {
    return this->fromConnectivity;
};
void DofNumbering::SetFromConnectivity(const bool& val){
    this->fromConnectivity = val;
};

void DofNumbering::ComputeNumberingFromConnectivity(UnstructuredMesh& mesh){
    this->size = mesh.GetNumberOfNodes();
    this->fromConnectivity = true;
    this->numbering.clear();
    for (auto& x : mesh.elements.storage){
        this->numbering[x.first] = x.second.GetConnectivityMatrix();
    }
    this->almanac.clear();
    for(CBasicIndexType i =0; i<this->size; i++ ){
        this->almanac[GetPointDofKey(i)] = i ;
    }
}

void DofNumbering::ComputeNumberingGeneral(UnstructuredMesh& mesh, Space& space, ElementFilterBase& elementFilter){
    CBasicIndexType& size = this->size;

    CBasicIndexType useddim = 0;
    CBasicIndexType cctt = 0;
    for (auto& x : mesh.elements.storage){
        MatrixID1 ids = elementFilter.GetIdsToTreat(mesh, x.first);
        CBasicIndexType  nel = static_cast<CBasicIndexType>(ids.rows());
        if (nel == 0) continue;
        const CBasicIndexType dim = ElementNames[x.first].dimension();
        useddim = std::max(useddim, dim);
        cctt += nel;
        const CBasicIndexType& nbOfShapeFuntions = space.GetNumberOfShapeFunctionsFor(x.first) ;
        const CBasicIndexType& nbOfElement = x.second.GetNumberOfElements();

        if ( !this->HasNumberingFor(x.first) ){
            this->InitNumberingFor(x.first, nbOfElement,nbOfShapeFuntions);
        }

        auto& localNumbering = this->GetNumberingFor(x.first);
        const auto& localspace = space.GetSpaceFor(x.first);

        const auto& localConnectivity = x.second.GetConnectivityMatrix();
        CBasicIndexType d;
        for(int j = 0 ; j < nbOfShapeFuntions; ++j){
            for(CBasicIndexType el=0;el < nel; ++el){

                const CBasicIndexType& currentElementId = ids(el,0);
                const auto& elcoon = localConnectivity.row(currentElementId);
                const auto& key = this->GetKeyFor(x.first,currentElementId, j,x.second,localspace.GetDofAttachment(j), elcoon );

                if (this->almanac.count(key)){
                    d = this->almanac[key];
                } else {
                    d = size;
                    this->almanac[key] = d;
                    ++size;
                }
                localNumbering(currentElementId,j) = d;
            };
        };
    }
    for (auto& x : mesh.elements.storage){
        if ( useddim <= ElementNames[x.first].dimension()) continue;

        MatrixID1 ids = elementFilter.GetIdsToTreatComplementaty(mesh, x.first);

        CBasicIndexType  nel = static_cast<CBasicIndexType>(ids.rows());
        if (nel == 0) continue;

        const auto& localspace = space.GetSpaceFor(x.first);
        const CBasicIndexType& nbOfShapeFuntions = space.GetNumberOfShapeFunctionsFor(x.first) ;
        const CBasicIndexType& nbOfElement = x.second.GetNumberOfElements();

        if ( !this->HasNumberingFor(x.first) ){
            this->InitNumberingFor(x.first, nbOfElement,nbOfShapeFuntions);
        }
        auto& localNumbering = this->GetNumberingFor(x.first);

        const auto localConnectivity = x.second.GetConnectivityMatrix();
        for(CBasicIndexType el=0;el < nel; ++el){
            const CBasicIndexType& currentElementId = ids(el,0);
            const auto& elcoon = localConnectivity.row(currentElementId);
            for(int j = 0 ; j < nbOfShapeFuntions; ++j){
                const auto& key = this->GetKeyFor(x.first,currentElementId, j,x.second,localspace.GetDofAttachment(j), elcoon );
                //CBasicIndexType d =-1;
                if (this->almanac.count(key)){
                    localNumbering(currentElementId,j) = this->almanac[key];
                }
            };
        };
    }
}

bool DofNumbering::HasNumberingFor(const std::string & elemtype){
    if (this->numbering.count(elemtype)>0){
        return this->numbering[elemtype].rows() >0 && this->numbering[elemtype].cols() >0;
    }
    return false;
}
void DofNumbering::InitNumberingFor(const std::string & elemtype, const CBasicIndexType& nbOfElement, const CBasicIndexType& nbOfShapeFuntions ){
    if (nbOfElement == 0){
        this->numbering.erase(elemtype);
        return ;
    }
    this->numbering[elemtype].resize(nbOfElement,nbOfShapeFuntions);
    this->numbering[elemtype].fill(-1);
}

MatrixIDD& DofNumbering::GetNumberingFor(const std::string & elemtype){
    if(this->HasNumberingFor(elemtype)){
        return this->numbering[elemtype];
    }
    std::cout << "element " << elemtype << std::endl;
    PRINTDEBUG(this->numbering[elemtype]);
    PRINTDEBUG("-*-*-*-*");
    std::cout << std::flush;
    throw "no numbering for" + elemtype;
}

CBasicIndexType DofNumbering::GetDofOfPoint(const CBasicIndexType& pid){
    if(this->almanac.count(GetPointDofKey(pid)))
        return this->almanac[GetPointDofKey(pid)];
    else
        return -1;
}

MatrixID1& DofNumbering::GetdoftopointLeft(){
    this->computeDofToPoint();
    return this->doftopointLeft;
}

MatrixID1& DofNumbering::GetdoftopointRight(){
    this->computeDofToPoint();
    return this->doftopointRight;
}
//
void DofNumbering::computeDofToPoint(){
    if (this->dofToPointComputed) return ;

    CBasicIndexType cpt =0;
    for(const  auto& dof : this->almanac){
        if (dof.first.IsPointDof()){
            cpt += 1;
        }
    }
    if (cpt == 0) return;
    this->doftopointLeft.resize(cpt,1);
    this->doftopointRight.resize(cpt,1);

    cpt =0;
    for(const  auto& dof : this->almanac){
        if (dof.first.IsPointDof() ){
            this->doftopointLeft(cpt,0) = dof.first.GetPoit();
            this->doftopointRight(cpt,0) = dof.second;
            cpt += 1;
        }
    }
    this->dofToPointComputed = true;
}

CBasicIndexType DofNumbering::GetDofOfkey(const DofKey& key){
    if(this->almanac.count(key))
        return this->almanac[key];
    else
        return -1;
}

CBasicIndexType DofNumbering::GetSizeOfDofToPoint()   {
    this->computeDofToPoint();
    return static_cast<CBasicIndexType>(this->doftopointLeft.rows());
}

CBasicIndexType DofNumbering::GetSizeOfDofToCell()  {
    return static_cast<CBasicIndexType>(this->doftocellLeft.rows());
}

void DofNumbering::computeDofToCell(UnstructuredMesh& mesh){
    if (this->dofToCellComputed) return ;

    CBasicIndexType dofcpt = mesh.elements.GetNumberOfElements();
    this->doftocellLeft.resize(dofcpt,1);
    this->doftocellRight.resize(dofcpt,1);
    dofcpt = 0;
    CBasicIndexType elemcpt = 0;
    for (auto& x : mesh.elements.storage){
        const auto& connectivityMatrix = x.second.GetConnectivityMatrix();
        const CBasicIndexType nbel = x.second.GetNumberOfElements();
        for(CBasicIndexType el = 0; el < nbel; ++el ){
            DofKey key = DofKey(x.first,connectivityMatrix.row(el),0);
            key.SortId0();
            CBasicIndexType dofnb = GetDofOfkey(key);
            if(dofnb == -1) continue;
            this->doftocellLeft(dofcpt,0) = el + elemcpt;
            this->doftocellRight(dofcpt,0) = dofnb;
            ++dofcpt;
        }
        elemcpt += nbel;
    }

    this->doftocellLeft.conservativeResize(dofcpt);//,Eigen::NoChange_t);
    this->doftocellRight.conservativeResize(dofcpt);//,Eigen::NoChange_t);
    this->dofToCellComputed = true;
}
//
MatrixID1& DofNumbering::GetdoftocellLeft(){
     if (!this->dofToCellComputed) throw "Call computeDofToCell first";
    return this->doftocellLeft;
}
//
MatrixID1& DofNumbering::GetdoftocellRight(){
    if (!this->dofToCellComputed) throw "Call computeDofToCell first";
    return this->doftocellRight;
}
//
std::string DofNumbering::ToStr(){
    std::string res = "";
    return res;

}

} // namespace BasicTools
