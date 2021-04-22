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

DofKey::DofKey(const std::string s,const int& id0, const int& id1): s(s), id1(id1){
    this->id0.resize(1,1);
    this->id0(0,0)= id0;
};

DofKey::DofKey(const std::string s,const INT_TYPE& id0, const int& id1): s(s), id1(id1){
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
               const int& m = this->id0.rows();
               for(int i =0;  i < m ; ++i){
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
INT_TYPE DofKey::GetPoit()const {
    return this->id0(0,0);
}
//------------------------------------------------------------
DofKey GetPointDofKey(const INT_TYPE&pid ) {
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
INT_TYPE DofNumbering::GetSize() const {
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
    for(INT_TYPE i =0; i<this->size; i++ ){
        this->almanac[GetPointDofKey(i)] = i ;
    }
}

void DofNumbering::ComputeNumberingGeneral(UnstructuredMesh& mesh, Space& space, ElementFilterBase& elementFilter){
    INT_TYPE& size = this->size;
    
    INT_TYPE useddim = 0;
    INT_TYPE cctt = 0;
    for (auto& x : mesh.elements.storage){
        MatrixID1 ids = elementFilter.GetIdsToTreat(mesh, x.first);
        INT_TYPE  nel = ids.rows();
        if (nel == 0) continue;
        const INT_TYPE dim = ElementNames[x.first].dimension();
        useddim = std::max(useddim, dim);
        cctt += nel;
        const int& nbOfShapeFuntions = space.GetNumberOfShapeFunctionsFor(x.first) ;
        const int& nbOfElement = x.second.GetNumberOfElements();
        
        if ( !this->HasNumberingFor(x.first) ){
            this->InitNumberingFor(x.first, nbOfElement,nbOfShapeFuntions);
        }
        
        auto& localNumbering = this->GetNumberingFor(x.first);
        const auto& localspace = space.GetSpaceFor(x.first);
        
        
        const auto& localConnectivity = x.second.GetConnectivityMatrix();
        INT_TYPE d;
        for(int j = 0 ; j < nbOfShapeFuntions; ++j){
            for(INT_TYPE el=0;el < nel; ++el){
                
                const INT_TYPE& currentElementId = ids(el,0);
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
        
        INT_TYPE  nel = ids.rows();
        if (nel == 0) continue;
        
        const auto& localspace = space.GetSpaceFor(x.first);
        const int& nbOfShapeFuntions = space.GetNumberOfShapeFunctionsFor(x.first) ;
        const int& nbOfElement = x.second.GetNumberOfElements();
        
        if ( !this->HasNumberingFor(x.first) ){
            this->InitNumberingFor(x.first, nbOfElement,nbOfShapeFuntions);
        }
        auto& localNumbering = this->GetNumberingFor(x.first);
        
        const auto localConnectivity = x.second.GetConnectivityMatrix();
        for(INT_TYPE el=0;el < nel; ++el){
            const INT_TYPE& currentElementId = ids(el,0);
            const auto& elcoon = localConnectivity.row(currentElementId);
            for(int j = 0 ; j < nbOfShapeFuntions; ++j){
                const auto& key = this->GetKeyFor(x.first,currentElementId, j,x.second,localspace.GetDofAttachment(j), elcoon );
                //INT_TYPE d =-1;
                if (this->almanac.count(key)){
                    localNumbering(currentElementId,j) = this->almanac[key];
                }
            };
        };
    }
}

bool DofNumbering::HasNumberingFor(const std::string & elemtype){
    return this->numbering.count(elemtype)>0;
}
void DofNumbering::InitNumberingFor(const std::string & elemtype, const long int& nbOfElement, const long int& nbOfShapeFuntions ){
    this->numbering[elemtype].resize(nbOfElement,nbOfShapeFuntions);
    this->numbering[elemtype].fill(-1);
}

MatrixIDD& DofNumbering::GetNumberingFor(const std::string & elemtype){
    if(this->numbering.count(elemtype)){
        return this->numbering[elemtype];
    }
    std::cout << "element " << elemtype << std::endl;
    PRINTDEBUG(this->numbering[elemtype]);
    PRINTDEBUG("-*-*-*-*");
    std::cout << std::flush;
    throw "no numbering for" + elemtype;
}

INT_TYPE DofNumbering::GetDofOfPoint(const INT_TYPE& pid){
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
    
    INT_TYPE cpt =0;
    for(const  auto dof : this->almanac){
        if (dof.first.IsPointDof()){
            cpt += 1;
        }
    }
    if (cpt == 0) return;
    this->doftopointLeft.resize(cpt,1);
    this->doftopointRight.resize(cpt,1);
    
    cpt =0;
    for(const  auto dof : this->almanac){
        if (dof.first.IsPointDof() ){
            this->doftopointLeft(cpt,0) = dof.first.GetPoit();
            this->doftopointRight(cpt,0) = dof.second;
            cpt += 1;
        }
    }
    this->dofToPointComputed = true;
}

INT_TYPE DofNumbering::GetDofOfkey(const DofKey& key){
    if(this->almanac.count(key))
        return this->almanac[key];
    else
        return -1;
}

INT_TYPE DofNumbering::GetSizeOfDofToPoint()   {
    this->computeDofToPoint();
    return this->doftopointLeft.rows();
}

INT_TYPE DofNumbering::GetSizeOfDofToCell()  {
    return this->doftocellLeft.rows();
}

void DofNumbering::computeDofToCell(UnstructuredMesh& mesh){
    if (this->dofToCellComputed) return ;
   
    INT_TYPE dofcpt = mesh.elements.GetNumberOfElements();
    this->doftocellLeft.resize(dofcpt,1);
    this->doftocellRight.resize(dofcpt,1);
    dofcpt = 0;
    INT_TYPE elemcpt = 0;
    //for (auto&& [elemtype, elemdata] : mesh.elements.storage){
    for (auto& x : mesh.elements.storage){
        const auto& connectivityMatrix = x.second.GetConnectivityMatrix();
        const INT_TYPE nbel = x.second.GetNumberOfElements();
        for(INT_TYPE el = 0; el < nbel; ++el ){
            DofKey key = DofKey(x.first,connectivityMatrix.row(el),0);
            key.SortId0();
            INT_TYPE dofnb = GetDofOfkey(key);
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
