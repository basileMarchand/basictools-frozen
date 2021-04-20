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
#include <ImplicitGeometry/ImplicitGeometryBase.h>

namespace BasicTools
{
    
class ElementFilterBase {
public:
    virtual const MatrixID1 GetIdsToTreat(UnstructuredMesh& mesh, const std::string& elemtype)  =0;
    const MatrixID1 GetIdsToTreatComplementaty(UnstructuredMesh& mesh, const std::string& elemtype);
};

class ElementFilterEvaluated : public ElementFilterBase {
    std::map<std::string,MatrixID1 > ids;
    std::map<std::string,INT_TYPE> nbelements;
public:
    const MatrixID1 GetIdsToTreat(UnstructuredMesh& mesh, const std::string& elemtype);
    //
    void SetIdsToTreatFor(const std::string& elemtype, const Eigen::Ref<const MatrixID1>& ids);
    //
    void Clear();
    //
    virtual std::string ToStr() const;
};

class ElementFilter : public ElementFilterBase {
public:
    enum ZONE {CENTER=0, ALLNODES, LEASTONENODE};
private:
    int dimensionality;
    ZONE zoneTreatment;
    std::vector<std::string> tags;
    std::vector<std::string> elementTypes;
    std::vector<std::shared_ptr<ImplicitGeometryBase> > zones;
    //
    bool CheckDimensionality(const std::string& elemtype,bool& active) const ;
    //
    bool CheckElementTypes(const std::string& elemtype,bool& active) const;
    //
    MatrixID1 CheckZones( UnstructuredMesh& mesh,  const std::string& elemtype, bool& active) const ;
public:
    //
    ElementFilter();
    //
    void AddTag(const std::string& tagName);
    void SetDimensionality(const int& dim);
    //
    virtual const  MatrixID1 CheckTags(Tags& tags,const INT_TYPE& ts, bool& active) const ;
    //
    virtual const  MatrixID1 GetIdsToTreat(UnstructuredMesh& mesh,  const std::string& elemtype) const ;
    //
    virtual std::string ToStr() const;
    void SetCENTER();
    void SetALLNODES();
    void SetLEASTONENODE();
};

// ----------------------------- Operator on ElementFilterBase
class ElementFilterIntersection : public ElementFilterBase {
    std::vector<ElementFilterBase*> storage;
public:
    const MatrixID1 GetIdsToTreat(UnstructuredMesh& mesh, const std::string& elemtype);
};


}// namespace BasicTools

