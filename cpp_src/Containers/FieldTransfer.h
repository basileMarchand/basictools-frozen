
#pragma once
#include <string>
#include <Containers\UnstructuredMesh.h>
#include <FE\DofNumbering.h>
#include <Containers\ElementFilter.h>


#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
// to store queries results
#include <vector>
// just for output
#include <iostream>
#include <boost/foreach.hpp>
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
typedef bg::model::point<BasicTools::CBasicFloatType, 3, bg::cs::cartesian> point;
//CBasicFloatType
//typedef bg::model::box<point> box;
//typedef std::pair<box, BasicTools::CBasicIndexType> value;
typedef std::pair<point, BasicTools::CBasicIndexType> value;






namespace BasicTools {


template<typename T>
double normsquared(const T& v){
    return v.dot(v);
};


enum TransferMethods{
    Nearest=0,
    Interp,
    Extrap,
    Clamp,
    ZeroFill
    };

class TransferClass {

public:
    TransferClass();
    std::string ToStr();
    void SetVerbose(bool verbose= false);
//    void SetSourceFEField();
    void SetSourceMesh(UnstructuredMesh* sourceMesh);
    void SetSourceSpace(const std::string& space);
    void SetSourceNumbering(DofNumbering* numbering);

    void SetTransferMethod(const std::string& method);
    std::string GetTransferMethod();
//
    void SetElementFilter(const ElementFilterEvaluated& filter);
    MAPSETGET_MatrixDDD(TargetPoints,targetPoints)

    void Compute();
//  GetProjectorOperator()
    std::vector<CBasicIndexType> rows;
    std::vector<CBasicIndexType> cols;
    std::vector<CBasicFloatType> data;
    CBasicIndexType nb_source_Dofs;
    CBasicIndexType nb_targetPoints;
private:
    bool verbose;
    int insideMethod;
    int outsideMethod;
    UnstructuredMesh* sourceMesh;
    Space sourceSpace;
    DofNumbering* sourceNumbering;
    ElementFilterEvaluated elementfilter;
    bool elementFilterSet;
    std::shared_ptr<MapMatrixDDD > targetPoints;

    bgi::rtree< value, bgi::quadratic<16> > nodeRTree;
    bgi::rtree< value, bgi::quadratic<16> > centerRTree;
    MatrixIDD dualGraph;
    MatrixID1 usedPoints;
    MatrixDDD cellsCenters;
    void ComputeBarycentricCoordinateOnElement(){}


};



template<typename A, typename B, typename C, typename D, typename S2, typename S3>
void  ComputeBarycentricCoordinateOnElementTemplate(const A& coordAtDofs, const B& localspace, const C& targetPoint, const std::string& elementType,
                                            bool& inside,
                                            S2& xietaphi,
                                            S3& xichietaClamped){


const int elemDim = localspace[this->geoSpaceNumber].dimensionality;
const int spaceDim = static_cast<int>(this->nodes->cols());



/*
    xietaphi = np.array([0.5]*spacedim)
    N = localspace.GetShapeFunc(xietaphi)
    currentPoint = N.dot(coordAtDofs)
    f = targetPoint - currentPoint

    for x in range(10):
        dN = localspace.GetShapeFuncDer(xietaphi)
        df_num = df(f,dN,coordAtDofs)
        H = ddf(f, xietaphi, dN, localspace.GetShapeFuncDerDer, coordAtDofs, linear)
        if spacedim == 2:
            dxietaphi = inv22(H).dot(df_num)
        elif spacedim == 3:
            dxietaphi = hdinv(H).dot(df_num)
        else:
            dxietaphi = df_num/H[0,0]
        xietaphi -= dxietaphi

        # if the cell is linear only one iteration is needed
        if linear :
            break

        N = localspace.GetShapeFunc(xietaphi)
        f = targetPoint - N.dot(coordAtDofs)

        if normsquared(dxietaphi) < 1e-3 and normsquared(f) < 1e-3 :
            break
    else:
        return None, xietaphi,localspace.ClampParamCoorninates(xietaphi)

    xichietaClamped = localspace.ClampParamCoorninates(xietaphi)
    # we treat very closes point as inside
    inside = normsquared(xichietaClamped-xietaphi) < 1e-5
    return inside, xietaphi, xichietaClamped
*/

}

}
