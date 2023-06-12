
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
    MatrixID1 status;

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
    //void ComputeBarycentricCoordinateOnElement(){}
    std::pair<ElementsContainer,int> GetElement(int enb );


};



}
