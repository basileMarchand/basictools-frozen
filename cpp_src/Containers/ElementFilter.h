//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//
#pragma once

#include <string>
#include <ostream>
#include <vector>
#include <memory>

#include <LinAlg/EigenTypes.h>
#include <LinAlg/EigenTools.h>
#include <Containers/ElementNames.h>
#include <Containers/UnstructuredMesh.h>
#include <Containers/UnstructuredMeshTools.h>


struct SpaceFunction{
    virtual MatrixDD1 Eval(const MatrixDDD& pos) const = 0;
};

class ElementFilterBase{
public:
    virtual const MatrixID1 GetIdsToTreat(NativeUnstructuredMesh& mesh, const std::string& elemtype)  =0;
    const MatrixID1 GetIdsToTreatComplementaty( NativeUnstructuredMesh& mesh, const std::string& elemtype);/*  {
        const MatrixID1& ids =  this->GetIdsToTreat(mesh, elemtype);
        //std::cout << "line:" << __LINE__  << ": ids "<< ids << std::endl;
        MatrixID1 cids;
        const long int nbelements = mesh.elements[elemtype].GetNumberOfElements();

        INT_TYPE cpt = 0;
        if(ids.rows()==0 ){
            PRINT(elemtype)
            PRINT(nbelements)
            cids.resize(nbelements-ids.size(),1);
            cids.setLinSpaced(nbelements,0,nbelements-1);
            return cids;
        };
        cids.resize(nbelements-ids.size(),1);

        PRINT( ": id "<< ids.rows() << "  cid" << cids.rows() <<  " " << cids.cols())

        for(INT_TYPE j = 0;j < ids(0,0); ++j ){
           cids(cpt,0) = j;
           //std::cout << __LINE__  << std::endl;
           ++cpt ;
        }

        for(INT_TYPE i = 0;i < ids.rows()-1; ++i ){
           for(INT_TYPE j = ids(i)+1;j < ids(i+1,0); ++j ){
               cids(cpt,0) = j;
               //std::cout << __LINE__  << std::endl;
               ++cpt ;
           }
        }

        for(INT_TYPE j = ids(ids.rows()-1,0)+1;j < nbelements; ++j ){
           cids(cpt,0) = j;
           //std::cout << __LINE__  << " " << cpt <<  " " << j << std::endl;
           ++cpt ;
        }

        return cids;
    }*/
    //virtual std::string ToStr() const = 0;
};

class ElementFilterIntersection : public ElementFilterBase {
    std::vector<ElementFilterBase*> storage;
public:
    const MatrixID1 GetIdsToTreat(NativeUnstructuredMesh& mesh, const std::string& elemtype); /* {
        MatrixID1 res = storage[0]->GetIdsToTreat(mesh, elemtype);
        for (unsigned int i = 0; i< this->storage.size(); ++i){
            res = Intersect1D(res,storage[i]->GetIdsToTreat(mesh, elemtype));
        }
        return res;
    };*/
};


class ElementFilterEvaluated : public ElementFilterBase {
    std::map<std::string,MatrixID1 > ids;
    std::map<std::string,INT_TYPE> nbelements;
public:
    const MatrixID1 GetIdsToTreat(NativeUnstructuredMesh& mesh, const std::string& elemtype)  {
        if( this->ids.count(elemtype) == 0 ) {
            this->ids[elemtype] = MatrixID1();//(nullptr, 0, 1);
            this->ids[elemtype].resize(0,1);
        }
        return this->ids[elemtype];
    }
    //
    template<typename T>
    void SetIdsToTreatFor(const std::string& elemtype, T& ids){
        // auto val = new std::shared_ptr<MapMatrixID1 >::element_type(ids.data(), ids.rows(), ids.cols());
        this->ids[elemtype] = ids;
    }
    //
    void Clear(){ this->ids.clear(); }
    //
    virtual std::string ToStr() const {
      std::string res =  "  ElementFilterEvaluated";
      return res;
    };
};


class ElementFilter : public ElementFilterBase {
public:
    enum ZONE{CENTER=0, ALLNODES, LEASTONENODE};
private:
    int dimensionality;
    ZONE zoneTreatment;
    std::vector<std::string> tags;
    std::vector<std::string> elementTypes;
    std::vector<std::shared_ptr<SpaceFunction> > zones;
    //
    bool CheckDimensionality(const std::string& elemtype,bool& active) const {
        if (this->dimensionality == -100){
            active = false;
            return true;
        }
        active = true;
        const ElementInfo& element = ElementNames[elemtype];
        const int eldim = element.dimension();
        if(this->dimensionality  >= 0){
            if (eldim != this->dimensionality){
                return false;
            }else{
                return true;
            }
        } else {
            if (eldim == -this->dimensionality){
                return false;
            } else {
                return true;
            }
        }
    };
    //
    bool CheckElementTypes(const std::string& elemtype,bool& active) const {
        if (this->elementTypes.size() == 0){
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
    MatrixID1 CheckZones( NativeUnstructuredMesh& mesh,  const std::string& elemtype, bool& active) const {

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
        for(unsigned int i = 0; i < this->zones.size(); ++i){
            MatrixBD1 res2;
            const SpaceFunction* zone = this->zones[i].get();
            if (this->zoneTreatment == CENTER) {
                res2 = (zone->Eval(centers).array()<=0);
            } else if (this->zoneTreatment == ALLNODES) {
                MatrixID1 z = (zone->Eval(mesh.GetNodesMatrix()).array()<=0).cast<INT_TYPE>().matrix();
                MatrixID1 AA = indexingm(z, elements.GetConnectivityMatrix()).rowwise().sum().eval();

                res2 =  AA.array() == (long int)ElementNames[elemtype].numberOfNodes;
            } else if (this->zoneTreatment == LEASTONENODE ) {
                MatrixID1 z = (zone->Eval(mesh.GetNodesMatrix()).array() <=0).cast<INT_TYPE>().matrix();
                res2 = indexingm(z,elements.GetConnectivityMatrix()).rowwise().sum().array() > 0;
            }
            res = (res && res2).eval();
        }
        return where(res);
    }
public:
    //
    ElementFilter(): dimensionality(-100), zoneTreatment(CENTER) {

    };
    //
    void AddTag(const std::string& tagName){
        this->tags.push_back(tagName);
    }
    //
    void SetDimensionality(const int& dim){
        if ((dim == -100 or dim > -4) and dim < 4){
            this->dimensionality = dim;
        }
        throw "Value of dim in SetDimensionality not valid: ";
    };
    //
    virtual const  MatrixID1 CheckTags(Tags& tags,const INT_TYPE& ts, bool& active) const {
        MatrixID1 res;

        if (this->tags.size() == 0){
            active = false;
            return res;
        }

        for(const std::string& tagname :  this->tags){
            if ( tags.Contains(tagname)) res = Union1D(res,tags[tagname].GetIdsMatrix());
        }
        return res;
    }
    //
    virtual const  MatrixID1 GetIdsToTreat(NativeUnstructuredMesh& mesh,  const std::string& elemtype) const {

        MatrixID1 res;
        bool active=false;
        bool bool_res = this->CheckDimensionality(elemtype,active);
        if (active and !bool_res) return res;

        bool_res = this->CheckElementTypes(elemtype,active);
        if (active and !bool_res) return res;


        res  = this->CheckTags(mesh.elements[elemtype].tags,mesh.elements[elemtype].GetNumberOfElements(),active);
        if (active and res.rows()==0) return res;

        bool activeZone = false;
        MatrixID1 res2 = CheckZones(mesh, elemtype, activeZone);
        if( active){
            if( activeZone){
                return Intersect1D(res,res2);
            } else {
                return res;
            }
        }else {
            if( activeZone){
                return res2;
            }else {
                res.setLinSpaced(mesh.elements[elemtype].GetNumberOfElements(), 0, mesh.elements[elemtype].GetNumberOfElements()-1);
                return res;
            }
        }
    }
    //
    virtual std::string ToStr() const{
        std::string res = "ElementFilter:" ;
        return res;
    };
};
