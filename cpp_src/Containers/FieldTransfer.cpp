
#include <Containers/FieldTransfer.h>
#include <FE/GeneratedSpaces.h>
#include <LinAlg/EigenTypes.h>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
// to store queries results
#include <vector>
// just for output
#include <iostream>
#include <boost/foreach.hpp>
/*
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
typedef bg::model::point<BasicTools::CBasicFloatType, 3, bg::cs::cartesian> point;
//CBasicFloatType
//typedef bg::model::box<point> box;
//typedef std::pair<box, BasicTools::CBasicIndexType> value;
typedef std::pair<point, BasicTools::CBasicIndexType> value;
*/

#include <array>
#include <iostream>
#include <string>

namespace BasicTools {


// to  https://www.boost.org/doc/libs/1_74_0/libs/geometry/doc/html/geometry/spatial_indexes/creation_and_modification.html
template<typename T, typename I >
auto BuildRTree(const T& nodes, const I& originalId ){

  //bgi::rtree< Value, index::linear<16> > rtree;
  bgi::rtree< value, bgi::quadratic<16> > rtree;
  //bgi::rtree< Value, index::rstar<16> > rtree;

  const auto nbnodes = nodes.rows();
  for(int n=0; n < nbnodes; ++n){
    rtree.insert(std::make_pair(point(nodes(n,0),nodes(n,1),nodes(n,2)), originalId[n]));
  }
  return rtree;
}

template<typename P>
BasicTools::CBasicIndexType FindNearest(bgi::rtree< value, bgi::quadratic<16> >& rtree, const P& pointCoord ){
  //std::vector<value> result_n;
  //rtree.query(bgi::nearest(point(pointCoord(0,0),pointCoord(0,1),pointCoord(0,2)), 1), std::back_inserter(result_n));
  //return result_n.front().second;
  //auto rtree.qbegin(bgi::nearest(point(pointCoord(0,0),pointCoord(0,1),pointCoord(0,2)), 1))
  std::cout << "Start FindNearest" << std::endl;
  std::cout << "pointCoord "  << pointCoord << std::endl;
  auto pp = point(pointCoord(0),pointCoord(1),pointCoord(2));
  std::cout << "second part" << std::endl;
  auto res  = rtree.qbegin(bgi::nearest(pp, 1));
  std::cout << "--- " <<   std::endl;
  rtree.qend();
  std::cout << "---- " <<   std::endl;
  if (res == rtree.qend()) std::cout << "error!!!!!!!!!!!!!" << std::endl;
  std::cout << "---- second  " << res->second <<  std::endl;
  return res->second;
}



TransferClass::TransferClass()
    : verbose(false), insideMethod(Interp), outsideMethod(Clamp), elementFilterSet(false){


    };

std::string TransferClass::ToStr() { return "TransferClass"; };
//

void TransferClass::SetVerbose(bool verbose) { this->verbose = verbose; };

std::string TransferClass::GetTransferMethod(){
  std::array<std::string, 5> options{
      {"Nearest", "Interp", "Extrap", "Clamp", "ZeroFill"}};
    return options[this->insideMethod] + "/" + options[this->outsideMethod];
}

void TransferClass::SetTransferMethod(const std::string& method) {
  std::string::size_type idx = 0;
  idx = method.find('/', 0);
  std::array<std::string, 5> options{
      {"Nearest", "Interp", "Extrap", "Clamp", "ZeroFill"}};

  if (idx == std::string::npos) {
    std::cout << "method : (" << method << ") not accepted";
    exit(1);
  }
  int i;
  for (i = 0; i < 5; ++i) {
    if (method.substr(0, idx) == options[i]) break;
  }
  if (i == 5) {
    std::cout << "method : (" << method << ") not accepted";
    exit(1);
  }
  this->insideMethod = i;

  for (i = 0; i < 5; ++i) {
    if (method.substr(idx + 1, method.size()) == options[i]) break;
  }
  if (i == 5) {
    std::cout << "method : (" << method << ") not accepted";
    exit(1);
  }
  this->outsideMethod = i;
};

void TransferClass::SetSourceMesh(UnstructuredMesh* sourceMesh){
  std::cout << "start SetSourceMesh" << std::endl;
  this->sourceMesh = sourceMesh;
  this->nodeRTree = BuildRTree(this->sourceMesh->GetNodesMatrix(), this->sourceMesh->GetOriginalIdsMatrix());
  //    kdt = cKDTree(iNodes)

    //TODO Compute KDTree on inputPoints (sourceMesh.nodes)
  std::cout << "End SetSourceMesh" << std::endl;

  // ComputeNodeToElementConnectivity
  this->dualGraph.resize(sourceMesh->GetNumberOfNodes(),200);
  this->usedPoints.resize(sourceMesh->GetNumberOfNodes(),1);
  this->usedPoints.fill(0);

  int cpt =0;
  for (auto& x : this->sourceMesh->elements.storage){
    const CBasicIndexType nbElements = x.second.GetNumberOfElements();
    const auto connMatrix = x.second.GetConnectivityMatrix();
    for(int i=0; i < nbElements; ++i){
      const auto row = connMatrix.row(i);
      for( int j=0; j < row.size(); ++j){
        this->dualGraph(j,this->usedPoints(j,0)) =  cpt;
        this->usedPoints[j] += 1;
      }
    }
  }
  // Generate Cells Centers
  cellsCenters.resize(this->sourceMesh->elements.GetNumberOfElements(),3) ;
  cellsCenters.fill(0);
  cpt =0;
  const auto nodes =   this->sourceMesh->GetNodesMatrix();

  for (auto& x : this->sourceMesh->elements.storage){
    const CBasicIndexType nbElements = x.second.GetNumberOfElements();
    const auto connMatrix = x.second.GetConnectivityMatrix();
    for( int c =0; c <3; c++){
      cellsCenters.block(cpt,c,nbElements,1) += indexingi(nodes,connMatrix,c);
    }
    cellsCenters /= connMatrix.cols();

    cpt += nbElements;
  }



}



void TransferClass::SetSourceSpace(const std::string& spaceName){
    this->sourceSpace = GetFESpaceFor(spaceName);
};

void TransferClass::SetSourceNumbering(DofNumbering* numbering){
    this->sourceNumbering = numbering;
};

void TransferClass::SetElementFilter(const ElementFilterEvaluated& filter){
    this->elementfilter = filter;
    this->elementFilterSet = true;

};

void TransferClass::Compute(){
  std::cout << "Start Compute" << std::endl;

  CBasicIndexType nbTargetPoints = (CBasicIndexType) this->targetPoints->rows();
  MatrixID1 ids;
  ids.resize(nbTargetPoints,1);

  CBasicIndexType id;
  CBasicIndexType col;
  this->rows.resize(0);
  this->cols.resize(0);
  this->data.resize(0);

  if( (this->insideMethod == 0) && (this->outsideMethod == 0) ){
    for( int tp = 0; tp < nbTargetPoints; ++tp ){
      std::cout << "Inside point " << tp << std::endl;
      id = FindNearest(this->nodeRTree, this->targetPoints->row(tp) );


      if( this->sourceNumbering->GetFromConnectivity() ){
        col  = id;
      } else {
        col = this->sourceNumbering->GetDofOfPoint(id);
      }
      rows.push_back(tp);
      cols.push_back(col);
      data.push_back(1.);
    }
    this->nb_source_Dofs = this->sourceNumbering->GetSize();
    this->nb_targetPoints = nbTargetPoints;
    std::cout << "End Compute" << std::endl;
    return;
  }
  //we build de Dual Connectivity


};




template <typename A, typename B, typename C, typename D, typename E,
          typename F>
auto ddf(const A& f, const B& xiEtaPhi, const B& dN,
         const D& GetShapeFuncDerDer, const E& coordAtDofs, const F& linear) {
  const auto dNX = dN.dot(coordAtDofs) return dNX.dot(dNX.T)
}

template <typename A, typename B, typename C>
auto df(const A& f, const B& dN, const C& coordAtDofs) {
  return -f.dot(dN.dot(coordAtDofs).transpose);
}
/*
    dNX = dN.dot(coordAtDofs)
    res = -f.dot(dNX.T)
    return res
*/

};  // namespace BasicTools