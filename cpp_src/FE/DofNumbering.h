//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//
#pragma once


#include <memory>

#include <LinAlg/EigenTypes.h>
#include <Containers/UnstructuredMesh.h>
#include <Containers/ElementFilter.h>
#include <FE/Space.h>


namespace BasicTools
{
    
extern const std::string Point;
extern const std::string Global;
extern const std::string IntegrationPoint;

class DofKey{
    std::string s;
    //INT_TYPE id0;
    MatrixID1 id0;
    INT_TYPE id1;
public:
    DofKey(const std::string s,const int& id0, const int& id1);
    DofKey(const std::string s,const INT_TYPE& id0, const int& id1);
    template<typename T>
    DofKey(const std::string s,const T& id0, const int& id1): s(s), id1(id1){
        this->id0 = id0;
    };
    bool operator<( DofKey const& b) const ;
    void SortId0();
    void Print() const ;
    //
    bool IsPointDof()const;
    INT_TYPE GetPoit()const ;
};

DofKey GetPointDofKey(const INT_TYPE&pid ) ;

class DofNumbering{
    INT_TYPE size;
    bool fromConnectivity;
    bool discontinuous;
    std::map<std::string,MatrixIDD> numbering;
    std::map<DofKey,INT_TYPE > almanac;
    MatrixID1 doftopointLeft;
    MatrixID1 doftopointRight;
    bool dofToPointComputed;
    MatrixID1 doftocellLeft;
    MatrixID1 doftocellRight;
    bool dofToCellComputed;
public:
    DofNumbering();
    INT_TYPE GetSize() const ;
    bool GetFromConnectivity() const ;
    void SetFromConnectivity(const bool& val);
    void ComputeNumberingFromConnectivity(UnstructuredMesh& mesh);
    //
    void ComputeNumberingGeneral(UnstructuredMesh& mesh, Space& space, ElementFilterBase& elementFilter);
    //
    bool HasNumberingFor(const std::string & elemtype);
    //
    void InitNumberingFor(const std::string & elemtype, const long int& nbOfElement, const long int& nbOfShapeFuntions );
    //
    MatrixIDD& GetNumberingFor(const std::string & elemtype);
    //
    INT_TYPE GetDofOfkey(const DofKey& key);
    //
    INT_TYPE GetDofOfPoint(const INT_TYPE& pid);
    INT_TYPE GetSizeOfDofToPoint();
    MatrixID1& GetdoftopointLeft();
    MatrixID1& GetdoftopointRight();
    void computeDofToPoint();
    //
    void computeDofToCell(UnstructuredMesh& mesh);
    INT_TYPE GetSizeOfDofToCell();
    MatrixID1& GetdoftocellLeft();
    MatrixID1& GetdoftocellRight();
    //
    
    template< typename S>
    DofKey GetKeyFor(const std::string & elemtype,const INT_TYPE& elid,const int& sf, const ElementsContainer& data,const DofAttachment& da, const S & elcoon) const {
    //DofKey GetKeyFor(const std::string & elemtype,const INT_TYPE& elid,const int& sf, const ElementsContainer& data,const DofAttachment& da, const Eigen::Ref<const MatrixIDD> & elcoon) const ;        
    if(da.entity == 'P'){
        const INT_TYPE& pid = da.entityNumber;
        return GetPointDofKey(elcoon(0,pid));
    } else if (da.entity == 'C'){
        DofKey res = DofKey(elemtype,elcoon,da.entityNumber);
        res.SortId0();
        return  res;
    } else if (da.entity == 'F'){
        const auto& face = ElementNames[elemtype].faces[da.entityNumber];
        MatrixID1 face_connectivity(face.second.rows(),1);
        for(Eigen::Index i = 0 ; i < face.second.rows() ; ++i)
            face_connectivity(i,0) = elcoon(0,face.second(i,0));
        DofKey res = DofKey(face.first.name,face_connectivity ,0);
        res.SortId0();
        //res.Print();
        return res;
    } else if (da.entity == 'E'){
        const auto& face2 = ElementNames[elemtype].faces2[da.entityNumber];
        MatrixID1 face2_connectivity(face2.second.rows(),1);
        for(Eigen::Index i = 0 ; i < face2.second.rows() ; ++i){
            face2_connectivity(i,0) = elcoon(0,face2.second(i,0));
        }
        DofKey res = DofKey(face2.first.name,face2_connectivity ,0);
        res.SortId0();
        //res.Print();
        return res;
    } else if(da.entity == 'G'){
        return DofKey(Global,0,da.extraKey);
    } else if(da.entity == 'I') {
        //std::cout << elcoon << std::endl;
        return DofKey(IntegrationPoint,elcoon,sf);
    }
    throw  true;
    return GetPointDofKey(-1);
}
    std::string ToStr();
};

} // namespace BasicTools
