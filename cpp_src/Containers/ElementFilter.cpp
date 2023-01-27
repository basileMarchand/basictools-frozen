//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//

#include <Containers/ElementFilter.h>
//////////////// ElementFilterBase

namespace BasicTools
{

const MatrixID1 ElementFilterBase::GetIdsToTreatComplementary( UnstructuredMesh& mesh, const std::string& elementType)  {
    const MatrixID1 ids =  this->GetIdsToTreat(mesh, elementType);
    //std::cout << "line:" << __LINE__  << ": ids "<< ids << std::endl;
    MatrixID1 cids;
    const CBasicIndexType nbElements = mesh.elements[elementType].GetNumberOfElements();

    CBasicIndexType cpt = 0;
    if(ids.rows()==0 ) {
        PRINTDEBUG(elementType)
        PRINTDEBUG(nbElements)
        cids.resize(nbElements-ids.size(),1);
        cids.setLinSpaced(nbElements,0,nbElements-1);
        return cids;
    };
    cids.resize(nbElements-ids.size(),1);

    PRINTDEBUG( ": id "<< ids.rows() << "  cid" << cids.rows() <<  " " << cids.cols())

    for(CBasicIndexType j = 0; j < ids(0,0); ++j ) {
        cids(cpt,0) = j;
        ++cpt ;
    }

    for(CBasicIndexType i = 0; i < ids.rows()-1; ++i ) {
        for(CBasicIndexType j = ids(i)+1; j < ids(i+1,0); ++j ) {
            cids(cpt,0) = j;
            ++cpt ;
        }
    }

    for(CBasicIndexType j = ids(ids.rows()-1,0)+1; j < nbElements; ++j ) {
        cids(cpt,0) = j;
        ++cpt ;
    }

    return cids;
}
//////////////// ElementFilterIntersection

const MatrixID1 ElementFilterIntersection::GetIdsToTreat(UnstructuredMesh& mesh, const std::string& elemtype) const {
    MatrixID1 res = storage[0]->GetIdsToTreat(mesh, elemtype);
    for (unsigned int i = 0; i< this->storage.size(); ++i) {
        res = Intersect1D(res,storage[i]->GetIdsToTreat(mesh, elemtype));
    }
    return res;
};
////////////// ElementFilterEvaluated
const MatrixID1 ElementFilterEvaluated::GetIdsToTreat(UnstructuredMesh& mesh, const std::string& elemtype)  const {
    const auto it =  this->ids.find(elemtype);
    if( it == this->ids.end() ) {
        return MatrixID1();
    }
    return it->second;
}
//
void ElementFilterEvaluated::SetIdsToTreatFor(const std::string& elemtype, const Eigen::Ref<const MatrixID1>& ids) {
    this->ids[elemtype] = ids;
}
//
void ElementFilterEvaluated::Clear() {
    this->ids.clear();
}
//
std::string ElementFilterEvaluated::ToStr() const {
    std::string res =  "  ElementFilterEvaluated";
    return res;
};
//////////////// ElementFilter
ElementFilter::ElementFilter(): dimensionality(-100), zoneTreatment(CENTER) {};

void ElementFilter::SetCENTER() {
    this->zoneTreatment = ElementFilter::CENTER;
};
void ElementFilter::SetALLNODES() {
    this->zoneTreatment = ElementFilter::ALLNODES;
};
void ElementFilter::SetLEASTONENODE() {
    this->zoneTreatment = ElementFilter::LEASTONENODE;
};

bool ElementFilter::CheckDimensionality(const std::string& elemtype,bool& active) const {
    if (this->dimensionality == -100) {
        active = false;
        return true;
    }
    active = true;
    const ElementInfo& element = ElementNames[elemtype];
    const int eldim = element.dimension();
    if(this->dimensionality  >= 0) {
        if (eldim != this->dimensionality) {
            return false;
        } else {
            return true;
        }
    } else {
        if (eldim == -this->dimensionality) {
            return false;
        } else {
            return true;
        }
    }
};
bool ElementFilter::CheckElementTypes(const std::string& elemtype,bool& active) const {
    if (this->elementTypes.size() == 0) {
        active = false;
        return true;
    }
    active = true;
    if(std::find(this->elementTypes.begin(), this->elementTypes.end(), elemtype) != this->elementTypes.end()) {
        return true;
    } else {
        return false;
    }
}
//
MatrixID1 ElementFilter::CheckZones( UnstructuredMesh& mesh,  const std::string& elemtype, bool& active) const {

    const ElementsContainer& elements = mesh.elements[elemtype];
    if(this->zones.size() == 0) {
        active = false;
        return MatrixID1();
    }

    active = true;

    MatrixBD1 res(elements.GetNumberOfElements(),1);
    res.fill(0);
    MatrixDDD centers;

    if (this->zoneTreatment == CENTER) {
        centers = GetElementsCenters(mesh.GetNodesMatrix(),elements);
    }
    for(unsigned int i = 0; i < this->zones.size(); ++i) {
        MatrixBD1 res2;
        const ImplicitGeometryBase* zone = this->zones[i].get();
        if (this->zoneTreatment == CENTER) {
            res2 = (zone->GetDistanceToPoints(centers).array()<=0);
        } else if (this->zoneTreatment == ALLNODES) {
            MatrixID1 z = (zone->GetDistanceToPoints(mesh.GetNodesMatrix()).array()<=0).cast<CBasicIndexType>().matrix();
            MatrixID1 AA = indexingm(z, elements.GetConnectivityMatrix()).rowwise().sum().eval();

            res2 =  AA.array() == (long int)ElementNames[elemtype].numberOfNodes;
        } else if (this->zoneTreatment == LEASTONENODE ) {
            MatrixID1 z = (zone->GetDistanceToPoints(mesh.GetNodesMatrix()).array() <=0).cast<CBasicIndexType>().matrix();
            res2 = indexingm(z,elements.GetConnectivityMatrix()).rowwise().sum().array() > 0;
        }
        res = (res && res2).eval();
    }
    return NonZero(res);
}
void ElementFilter::AddTag(const std::string& tagName) {
    this->tags.push_back(tagName);
}
//
void ElementFilter::SetDimensionality(const int& dim) {
    if ((dim == -100 || dim > -4) && dim < 4) {
        this->dimensionality = dim;
    }
    throw "Value of dim in SetDimensionality not valid: ";
};
//
const  MatrixID1 ElementFilter::CheckTags(Tags& tags,const CBasicIndexType& ts, bool& active) const {
    MatrixID1 res;

    if (this->tags.size() == 0) {
        active = false;
        return res;
    }

    for(const std::string& tagname :  this->tags) {
        if ( tags.Contains(tagname)) res = Union1D(res,tags[tagname].GetIdsMatrix());
    }
    return res;
}
//
const  MatrixID1 ElementFilter::GetIdsToTreat(UnstructuredMesh& mesh,  const std::string& elemtype) const {

    MatrixID1 res;
    bool active=false;
    bool bool_res = this->CheckDimensionality(elemtype,active);
    if (active && !bool_res) return res;

    bool_res = this->CheckElementTypes(elemtype,active);
    if (active && !bool_res) return res;


    res  = this->CheckTags(mesh.elements[elemtype].tags,mesh.elements[elemtype].GetNumberOfElements(),active);
    if (active && res.rows()==0) return res;

    bool activeZone = false;
    MatrixID1 res2 = CheckZones(mesh, elemtype, activeZone);
    if( active) {
        if( activeZone) {
            return Intersect1D(res,res2);
        } else {
            return res;
        }
    } else {
        if( activeZone) {
            return res2;
        } else {
            res.setLinSpaced(mesh.elements[elemtype].GetNumberOfElements(), 0, mesh.elements[elemtype].GetNumberOfElements()-1);
            return res;
        }
    }
}
//
std::string ElementFilter::ToStr() const {
    std::string res = "ElementFilter:" ;
    return res;
};


}// namespace BasicTools
