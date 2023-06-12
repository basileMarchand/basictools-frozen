
#include <Containers/FieldTransfer.h>
#include <FE/GeneratedSpaces.h>
#include <LinAlg/EigenTypes.h>
#include <LinAlg/BasicOperations.h>
#include <Containers/ElementUtilities.h>

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

struct CEBCResult {
  MatrixD31 distv;
  MatrixD31 bary;
  MatrixD31 baryClamped;
  bool inside;
  bool error;
  CEBCResult(): inside(0) ,error(false){
    distv.fill(0);
    bary.fill(0);
    baryClamped.fill(0);
  };
};

template<typename C>
CEBCResult ComputeInterpolationExtrapolationsBarycentricCoordinates(const MatrixDDD& xcoor, const BasicTools::ElementSpace& localspace, const C& targetPoint);

// to  https://www.boost.org/doc/libs/1_74_0/libs/geometry/doc/html/geometry/spatial_indexes/creation_and_modification.html
template<typename T, typename I >
auto BuildRTree(const T& nodes, const I& originalId ){

  //bgi::rtree< Value, index::linear<16> > rtree;
  bgi::rtree< value, bgi::quadratic<16> > rtree;
  //bgi::rtree< Value, index::rstar<16> > rtree;

  const auto nbnodes = nodes.rows();
  for(int n=0; n < nbnodes; ++n){
    rtree.insert(std::make_pair(point(nodes(n,0),nodes(n,1),nodes(n,2)), n)); //originalId[n]
  }
  return rtree;
}

template<typename P>
BasicTools::CBasicIndexType FindNearest(bgi::rtree< value, bgi::quadratic<16> >& rtree, const P& pointCoord ){
  //std::vector<value> result_n;
  //rtree.query(bgi::nearest(point(pointCoord(0,0),pointCoord(0,1),pointCoord(0,2)), 1), std::back_inserter(result_n));
  //return result_n.front().second;
  //auto rtree.qbegin(bgi::nearest(point(pointCoord(0,0),pointCoord(0,1),pointCoord(0,2)), 1))
  //std::cout << "Start FindNearest" << std::endl;
  //std::cout << "pointCoord "  << pointCoord << std::endl;
  auto pp = point(pointCoord(0),pointCoord(1),pointCoord(2));
  //std::cout << "second part" << std::endl;
  auto res  = rtree.qbegin(bgi::nearest(pp, 1));
  //std::cout << "--- " <<   std::endl;
  rtree.qend();
  //std::cout << "---- " <<   std::endl;
  if (res == rtree.qend()) std::cout << "error!!!!!!!!!!!!!" << std::endl;
  //std::cout << "---- second  " << res->second <<  std::endl;
  return res->second;
}



TransferClass::TransferClass()
    : verbose(false), insideMethod(Interp), outsideMethod(Clamp), elementFilterSet(false){
      nb_source_Dofs = 0;
      nb_targetPoints = 0;
    };

std::string TransferClass::ToStr() { return "TransferClass"; };
//

void TransferClass::SetVerbose(bool verbose) { this->verbose = verbose; };
MatrixID1& TransferClass::GetStatus(){return this->status;}

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
 //std::cout << "-------------" << this->insideMethod << this->outsideMethod << std::endl;
};

void TransferClass::SetSourceMesh(UnstructuredMesh* sourceMesh){
 //std::cout << "start SetSourceMesh" << std::endl;
  this->sourceMesh = sourceMesh;
  this->nodeRTree = BuildRTree(this->sourceMesh->GetNodesMatrix(), this->sourceMesh->GetOriginalIdsMatrix());
  //    kdt = cKDTree(iNodes)

    //TODO Compute KDTree on inputPoints (sourceMesh.nodes)
 //std::cout << "End SetSourceMesh" << std::endl;

  // ComputeNodeToElementConnectivity
  this->dualGraph.resize(sourceMesh->GetNumberOfNodes(),200);
  this->dualGraph.fill(0);
  this->usedPoints.resize(sourceMesh->GetNumberOfNodes(),1);
  this->usedPoints.fill(0);

  int cpt =0;
  for (auto& x : this->sourceMesh->elements.storage){
    const CBasicIndexType nbElements = x.second.GetNumberOfElements();
    const auto connMatrix = x.second.GetConnectivityMatrix();
    for(int i=0; i < nbElements; ++i){
      const auto row = connMatrix.row(i);
      //std::cout << " row " << i << " " << row << std::endl;
      for( int j=0; j < row.size(); ++j){
        this->dualGraph(row(0,j),this->usedPoints(row(0,j),0)) =  cpt;
        this->usedPoints(row(0,j),0) += 1;
      }
      ++cpt;
    }
  }
  //std::cout << this->dualGraph << std::endl;
  //std::cout << this->usedPoints << std::endl;

  // Generate Cells Centers
  int nbElements = this->sourceMesh->elements.GetNumberOfElements();
  cellsCenters.resize(nbElements,3) ;
  cellsCenters.fill(0);
  cpt =0;
  const auto nodes =   this->sourceMesh->GetNodesMatrix();

  for (auto& x : this->sourceMesh->elements.storage){
    const CBasicIndexType nbElements = x.second.GetNumberOfElements();
    const auto connMatrix = x.second.GetConnectivityMatrix();
    for( int c =0; c <3; c++){
      for( int j=0; j < connMatrix.cols(); ++j){
        const MatrixDDD& res =  indexingi(nodes,connMatrix.col(j),c);
        cellsCenters.block(cpt,c,nbElements,1) += res;
      }
      cellsCenters.block(cpt,c,nbElements,1) *= 1./connMatrix.cols();
    }

    cpt += nbElements;
  }

  this->centerRTree = BuildRTree(cellsCenters,  MatrixID1::LinSpaced(nbElements,0,nbElements-1));


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

std::pair<ElementsContainer,int> TransferClass::GetElement(int enb){
  for (auto& x : this->sourceMesh->elements.storage){
    if ( enb < x.second.GetNumberOfElements()) {
        return std::make_pair(x.second, enb );
    } else {
       enb -= x.second.GetNumberOfElements();
    }
  }
  throw "Element not found";
}



struct BestCandidate{
  MatrixDD1 shapeFunc;
  MatrixDD1 shapeFuncClamped;
  MatrixID1 localnumbering;
  int memlenb;
  bool error;
  BestCandidate():error(false){};
};

void TransferClass::Compute(){
  //std::cout << "Start Compute**" << std::endl;

  CBasicIndexType nbTargetPoints = (CBasicIndexType) this->targetPoints->rows();

  this->nb_source_Dofs = this->sourceNumbering->GetSize();
  this->nb_targetPoints = nbTargetPoints;

  MatrixID1 ids;
  ids.resize(nbTargetPoints,1);

  CBasicIndexType id;
  CBasicIndexType col;
  // 30 to be sure to hold exa27 coefficients
  this->rows.reserve(nbTargetPoints*30);
  this->cols.reserve(nbTargetPoints*30);
  this->data.reserve(nbTargetPoints*30);

  this->rows.clear();
  this->cols.clear();
  this->data.clear();

  this->status.resize(nbTargetPoints,1);
  this->status.fill(0);

  const auto originalIdsPoints =  this->sourceMesh->GetOriginalIdsMatrix();
  const auto& sourceNodes = this->sourceMesh->GetNodesMatrix();

  if( (this->insideMethod == 0) && (this->outsideMethod == 0) ){
    for( int tp = 0; tp < nbTargetPoints; ++tp ){
      //std::cout << "Inside point " << tp << std::endl;
      id = FindNearest(this->nodeRTree, this->targetPoints->row(tp) );
      id = originalIdsPoints(id,0);

      if( this->sourceNumbering->GetFromConnectivity() ){
        col  = id;
      } else {
        col = this->sourceNumbering->GetDofOfPoint(id);
      }
      rows.push_back(tp);
      cols.push_back(col);
      data.push_back(1.);
    }

    //std::cout << "End Compute c++ nearest" << std::endl;
    return;
  }
  //we build de Dual Connectivity

  if (this->verbose){
    std::cout << "Starting c++ " << std::endl;
  }
  std::vector<CBasicIndexType> potentialElements;
  potentialElements.reserve(50);

  std::vector<CBasicFloatType> potentialElementsDistances;
  potentialElementsDistances.reserve(50);

  std::vector<CBasicIndexType> potentialElementIndex;
  potentialElementIndex.reserve(50);

  CBasicFloatType distmem;
  //CEBCResult resultMem;
  BestCandidate bestCandidate;
  for(int p= 0; p < nbTargetPoints; ++p){
    //if (this->verbose){
    //  std::cout << "|" ;
    //}
    const MatrixD1D TP = this->targetPoints->row(p);             // target point position
    //std::cout << " TP : "  << TP << std::endl;
    const auto cp_id = FindNearest(this->nodeRTree, TP );   // closest point id
    //const auto CP =  sourceNodes.row(cp_id);                // closest point position
    //std::cout << " cp_id : "  << cp_id << std::endl;

    const auto ce_id = FindNearest(this->centerRTree, TP );  // closest element id
    auto elementData_and_id = GetElement(ce_id);
    //std::cout << " ce_id : "  << ce_id << std::endl;
    //Element connected to the closest point
    potentialElements.clear();

    for(int dg_cpt=0; dg_cpt< this->usedPoints[cp_id] ; ++dg_cpt ){
      potentialElements.push_back(this->dualGraph(cp_id, dg_cpt) );
    }
    //Elements connected to the closest element (bases on the element center)
    for(int elempoint=0; elempoint < elementData_and_id.first.GetConnectivityMatrix().cols(); ++elempoint){
        const int lcoon = elementData_and_id.first.GetConnectivityMatrix()(elementData_and_id.second,elempoint);
        for(int dg_cpt=0; dg_cpt< this->usedPoints(lcoon,0) ; ++dg_cpt ){
          potentialElements.push_back(this->dualGraph(lcoon, dg_cpt) );
        }
    }

    std::vector<CBasicIndexType>::iterator it;
    std::sort(potentialElements.begin(), potentialElements.end());
    it = std::unique(potentialElements.begin(), potentialElements.end());
    potentialElements.erase(it,potentialElements.end());

    //std::cout << "potentialElementsD : " ;
    //for(long unsigned int i =0; i<potentialElements.size(); ++i){std::cout << potentialElements[i]; }
    //std::cout << std::endl; ;

    // compute distance**2 to elements
    // for the moment we use the distance to the center, this gives a good estimate
    // of the order to check the elements
    potentialElementsDistances.clear();
    potentialElementIndex.clear();
    for( long unsigned int potential_element_cpt=0; potential_element_cpt< potentialElements.size(); ++potential_element_cpt){
        const int e = potentialElements[potential_element_cpt];
        const auto diff = this->cellsCenters.row(e)-TP;
        potentialElementsDistances.push_back(diff.squaredNorm() );
        potentialElementIndex.push_back(potential_element_cpt);
    }
    // order the element to test, closest element first
    std::sort( potentialElementIndex.begin(),potentialElementIndex.end(), [&](int i,int j){
      return potentialElementsDistances[i]<potentialElementsDistances[j];
      });

    //std::cout << "potentialElementsDistances : " ;
    //for(long unsigned int i =0; i<potentialElementsDistances.size(); ++i){std::cout << potentialElementsDistances[i] << " "; } std::cout << std::endl; ;
    //for(long unsigned int i =0; i<potentialElementIndex.size(); ++i){std::cout << potentialElementIndex[i] << " "; } std::cout << std::endl; ;

    distmem = 1e10;
    //resultMem.rese()
    //computation for the real distance
    long unsigned int tecpt;
    for(tecpt=0; tecpt< potentialElements.size(); ++tecpt){
      const CBasicIndexType potential_element_cpt = potentialElementIndex[tecpt];
      const CBasicIndexType peGlobalId = potentialElements[potential_element_cpt];
      elementData_and_id = GetElement(peGlobalId);
      const int peLocalId = elementData_and_id.second;
      auto& localnumbering = this->sourceNumbering->GetNumberingFor(elementData_and_id.first.GetElementType());
      const ElementSpace& localspace = this->sourceSpace.GetSpaceFor(elementData_and_id.first.GetElementType());

      auto const pEConnectivity = elementData_and_id.first.GetConnectivityMatrix().row(peLocalId);
      const auto xcoor = sourceNodes(pEConnectivity,Eigen::all);
      //std::cout << " xcoor " << xcoor << std::endl;


      const CEBCResult result = ComputeInterpolationExtrapolationsBarycentricCoordinates(xcoor, localspace, TP);

      //update the distance**2 with a *exact* distance
      potentialElementsDistances[potential_element_cpt] = result.distv.squaredNorm();

     //std::cout << "------------   bary " << result.bary << std::endl;
     //std::cout << "bary Clamped " << result.baryClamped << std::endl;
     //std::cout << "inside " << result.inside << std::endl;
     //std::cout << "error " << result.error << std::endl;
     //std::cout << "real distance " << potentialElementsDistances[potential_element_cpt] << std::endl;
      //exit(1);
      //#compute shape function of the incomming space using the xi eta phi
      const auto shapeFunc = localspace.GetValOfShapeFunctionsAt(result.bary);
      const auto shapeFuncClamped =  localspace.GetValOfShapeFunctionsAt(result.baryClamped);
      // need to add a tolerance over the distv (real distance). this is needed because
      // we can have a point that the projection is inside an element (surface or line)
      // but not on the surface/line. Need to add a better test
      if (result.inside && (result.distv.squaredNorm() <= 1e-10) ){
        //std::cout << " ---------------------------------------- inside" << std::endl;
        bestCandidate.shapeFunc =  shapeFunc;
        bestCandidate.shapeFuncClamped =  shapeFuncClamped;
        bestCandidate.memlenb = elementData_and_id.first.GetOriginalIdsMatrix()(peLocalId,0);
        bestCandidate.localnumbering = localnumbering.row(bestCandidate.memlenb);
        status[p] = 1;
        break;
        }
      // store the best element (closest)
      if(potentialElementsDistances[potential_element_cpt] < distmem ){
        //std::cout << " not found yet -------------------------------------"  << std::endl;
        distmem = potentialElementsDistances[potential_element_cpt];
        bestCandidate.error = result.error;
        bestCandidate.shapeFunc =  shapeFunc;
        bestCandidate.shapeFuncClamped =  shapeFuncClamped;
        bestCandidate.memlenb = elementData_and_id.first.GetOriginalIdsMatrix()(peLocalId,0);
        bestCandidate.localnumbering = localnumbering.row(bestCandidate.memlenb);
      }
    }
    //# we are outside
   //std::cout << tecpt << " -*-*-*-* " <<  status[p] << std::endl;
    if( tecpt == potentialElements.size()){
      // no element found
      // or outsideMethod == Nearest
      // or extrapolation but error
     //std::cout << " ----------------------------------------------------------------------------------- "  << std::endl;
     //std::cout << " distmem "  << distmem << std::endl;
     //std::cout << " this->outsideMethod "  << this->outsideMethod << std::endl;
     //std::cout << " bestCandidate.error "  << bestCandidate.error << std::endl;
      //std::cout << " "  << << std::endl;

      if (distmem == 1e10 || this->outsideMethod == TransferMethods::Nearest || (outsideMethod == TransferMethods::Extrap && bestCandidate.error) ) {
       //std::cout << " AAAAAAAAA "  << std::endl;
        rows.push_back(p);
        cols.push_back(this->sourceNumbering->GetDofOfPoint(cp_id));
        data.push_back(1.);
        continue;
        }
      if (outsideMethod == TransferMethods::Extrap && bestCandidate.error == false){
       //std::cout << " BBBBBBBBB "  << std::endl;
        rows.insert(rows.end(),bestCandidate.shapeFunc.size(), p);
        cols.insert(cols.end(),bestCandidate.localnumbering.begin(), bestCandidate.localnumbering.end() );
        data.insert(data.end(), bestCandidate.shapeFunc.begin(), bestCandidate.shapeFunc.end()) ;
       //std::cout << " bestCandidate.shapeFunc "  << bestCandidate.shapeFunc << std::endl;

        status[p] = 2;
      } else if( outsideMethod == TransferMethods::Clamp){
       //std::cout << " CCCCCCCCCCCCCCCC "  << std::endl;

        rows.insert(rows.end(),bestCandidate.shapeFuncClamped.size(), p);
        cols.insert(cols.end(),bestCandidate.localnumbering.begin(), bestCandidate.localnumbering.end() );
        data.insert(data.end(), bestCandidate.shapeFuncClamped.begin(), bestCandidate.shapeFuncClamped.end()) ;
        status[p] = 3;
      } else if( outsideMethod == TransferMethods::ZeroFill){
       //std::cout << " DDDDDDDDDDDDDD "  << std::endl;

        //zero fill
        status[p] = 4;
      }
    } else {
      // found one element
     //std::cout << " adding point  : " <<  rows.size() << std::endl;
     //std::cout << " lenb------->" <<  bestCandidate.memlenb <<  "<-"<< std::endl;
     //std::cout << " col------->" <<  bestCandidate.localnumbering <<  "<-"<< std::endl;
     //std::cout << " adding point  : " <<  rows.size() << std::endl;
      rows.insert(rows.end(),bestCandidate.shapeFunc.size(), p);
      cols.insert(cols.end(),bestCandidate.localnumbering.begin(), bestCandidate.localnumbering.end() );
      data.insert(data.end(), bestCandidate.shapeFunc.begin(), bestCandidate.shapeFunc.end()) ;
    }

  }
 //std::cout << "  rows  : " <<  rows.size() << std::endl;
 //std::cout << "  cols  : " <<  cols.size() << std::endl;
 //std::cout << "  data  : " <<  data.size() << std::endl;
 //std::cout << "End Compute c++ rest" << std::endl;

};


template <typename B, typename C, typename D, typename E>
auto ddf(const MatrixD1D& f, const B& xiEtaPhi, const C& dN, const D& coordAtDofs, const E& linear) {
  const auto dNX = dN*(coordAtDofs) ;
  return dNX*(dNX.transpose());
}

template <typename B, typename C>
auto df(const MatrixD1D& f, const B& dN, const C& coordAtDofs) {
  //const auto temp = (dN*coordAtDofs).transpose();
 //std::cout << "temp" << temp << std::endl;
  return -f*((dN*coordAtDofs).transpose());
}

template<typename C>
CEBCResult ComputeInterpolationExtrapolationsBarycentricCoordinates(const MatrixDDD& xcoor, const BasicTools::ElementSpace& localspace, const C& targetPoint){

  if (localspace.dimensionality == 0){
    CEBCResult result;
    result.error = false;
    result.distv = (xcoor.row(0) - targetPoint).transpose();
    result.inside = result.distv.isZero();
    return result;
  }

  // we compute the baricentric coordinate of the target point (TP) on the current element
  CEBCResult temp_result = ComputeBarycentricCoordinateOnElement(xcoor, localspace, targetPoint);
  //std::cout << "---***---------" << std::endl;
  //std::cout << "inside " << temp_result.inside << std::endl;
  if (temp_result.inside){
    // the point is inside, we compute the distance vector
    temp_result.distv = (localspace.GetValOfShapeFunctionsAt(temp_result.bary).transpose()*xcoor).row(0)- targetPoint;
    return temp_result;
  } else {
    // compute distance to the corner points
    MatrixDD1 dist2 = (xcoor.rowwise() - targetPoint).array().square().rowwise().sum();
    //std::cout << " ======xcoor=====  " << xcoor << std::endl;
    //std::cout << " ======targetPoint=====  " << targetPoint << std::endl;
    //std::cout << " ===========  " << dist2 << std::endl;
    //std::cout << " =====-======  " << (xcoor.rowwise() - targetPoint) << std::endl;
    Eigen::Index minRow, minCol;
    dist2.minCoeff(&minRow, &minCol);
    MatrixID1 mask;
    mask.resize(dist2.rows(),1);
    mask.fill(0);
    mask(minRow,0) = 1;
    //std::cout << " =====mask=====  " << mask << std::endl;

    for(int faceNb =1 ; faceNb < 4 ; ++faceNb){
      //std::cout << " ===== !!!!! 1111 !!!!!!!!! =====  " << std::endl;
      const std::vector<std::pair<ElementInfo,MatrixID1> >& faces  = ElementNames[localspace.elementType].GetFacesLevel(faceNb);
      //std::cout << " ===== !!!!!!!!!!!!!! =====  " << std::endl;
      if (faces.size() == 0){
        //std::cout << " ===== Break =====  " << std::endl;

        break;
      }
      if ( faces[0].first.name == Point_1 ){
        // in this case we have in minRow the closes point
        CEBCResult result;
        result.error = false;
        result.bary = temp_result.bary;
        result.baryClamped.block(0,0,localspace.posN.cols(),1) =  localspace.posN.block(minRow,0,1,localspace.posN.cols()).transpose();// .row(minRow);
        result.distv =  (xcoor.row(minRow) - targetPoint).transpose();
        result.inside = false;
       //std::cout << " ===== point =====  " << std::endl;
        return result;
      }

      //faceInside, faceDistv, fbaryClamped ;
      CEBCResult result = ComputeInterpolationCoefficients(xcoor, localspace, targetPoint, mask,  faces );
      result.bary = temp_result.bary;
      if(result.inside){
        result.inside = false;
        return result;
      }
    }
    assert(0);
  }

  CEBCResult result;
  result.error = true;
  return result;
}

template<typename A, typename C >
CEBCResult ComputeInterpolationCoefficients(const A& xcoor,
                                            const BasicTools::ElementSpace& localspace,
                                            const C& targetPoint,
                                            const MatrixID1& mask,
                                            const std::vector<std::pair<ElementInfo,MatrixID1> >& faces ){

    CBasicFloatType faceDistMem = std::numeric_limits<CBasicFloatType>::max();
    CEBCResult result;
    result.error = true;
    result.inside = false;

    for( long unsigned int cpt =0; cpt < faces.size(); ++cpt){
      const std::string& faceElementType = faces[cpt].first.name;
      const MatrixID1& faceLocalConnectivity= faces[cpt].second;
      // work only on element touching the mask
      //if ( indexingi(mask,faceLocalConnectivity,0).isZero()){
      if ( mask(faceLocalConnectivity,Eigen::all).isZero()){

        continue;
      }
      const BasicTools::ElementSpace& lspace = GetFESpaceFor("LagrangeSpaceGeo").GetSpaceFor(faceElementType);
      const MatrixDDD localxcoor = xcoor(faceLocalConnectivity,Eigen::all);
      CEBCResult temp_result = ComputeBarycentricCoordinateOnElement(localxcoor, lspace, targetPoint);

      const MatrixDD1 faceDistv = ((lspace.GetValOfShapeFunctionsAt(temp_result.bary).transpose()*localxcoor).row(0) - targetPoint).transpose();

      CBasicFloatType  faceDist = faceDistv.squaredNorm();
      if (temp_result.inside && faceDist < faceDistMem){

          faceDistMem = faceDist;
          result.error = false;
          result.inside = true;
          result.distv = faceDistv;
          result.bary = temp_result.bary;
          //result.baryClamped =  temp_result.baryClamped;
          //flocalspace = posspace[faceElementTypeMem]
          const auto fshapeFunc = lspace.GetValOfShapeFunctionsAt(result.bary);
          const int elemDim = localspace.dimensionality;
          result.baryClamped.block(0,0,elemDim,1) = (fshapeFunc.transpose()*localspace.posN(faceLocalConnectivity,Eigen::all)).transpose();
      }
    }
    return result;
};

template<typename A, typename C>
CEBCResult ComputeBarycentricCoordinateOnElement(const A& xcoor, const BasicTools::ElementSpace& localspace, const C& targetPoint){
  /*compute element iternal coordinate on target point (targetPoint)*/
  CEBCResult result;
  const int elemDim = localspace.dimensionality;
  const int spaceDim = localspace.GetDimensionality();
  const int nbShapeFunctions = localspace.GetNumberOfShapeFunctions();
  //std::cout << "elemDim " << elemDim << std::endl;
  //std::cout << "spaceDim " << spaceDim << std::endl;
  //std::cout << "nbShapeFunctions " << nbShapeFunctions << std::endl;
  const bool linear = ElementNames[localspace.elementType].linear;
  //std::cout << "linear " << linear << std::endl;

  MatrixD13 xietaphi;// xietaphi.resize(1,3);
  xietaphi.fill(0.5);

  MatrixD1D N = localspace.GetValOfShapeFunctionsAt(xietaphi).transpose().row(0);
  MatrixD1D currentPoint = N*xcoor;
  MatrixD1D f =  currentPoint - targetPoint ;
  MatrixDD1 dxietaphi;
  int x=0;
  for(; x < 10; ++x){
    //std::cout << "----------- in iteration "<< x << std::endl;
    MatrixDDD dN = localspace.GetValOfShapeFunctionsDerAt(xietaphi).block(0,0,elemDim,nbShapeFunctions);
    //std::cout << " dN "<< MatrixDDD(dN) << std::endl;
    //std::cout << " f "<< MatrixDDD(f) << std::endl;
    //std::cout << " xcoor "<< MatrixDDD(xcoor) << std::endl;
    MatrixDDD df_num = df(-f,dN,xcoor).transpose();
    //std::cout << " df_num "<< df_num << std::endl;
    MatrixDDD H = ddf(-f, xcoor, dN,  xcoor, linear);
    //std::cout << " H "<< H << std::endl;
    if(spaceDim == 2){
      MatrixDDD Hinv;
      Hinv.resize(2,2);
      CBasicFloatType det;
      inv22(H,Hinv,det);
      dxietaphi = Hinv*df_num;
    } else if( spaceDim == 3 ){
      MatrixDDD Hinv;
      Hinv.resize(3,3);
      CBasicFloatType det;
      inv33(H,Hinv,det);
      dxietaphi = Hinv*df_num;
    } else{
      dxietaphi = df_num/H(0,0);
    }

    //std::cout << " dxietaphi "<< dxietaphi << std::endl;

    xietaphi.block(0,0,elemDim,1) -= dxietaphi.block(0,0,elemDim,1)    ;
    //# if the cell is linear only one iteration is needed
    if (linear){
      f *= 0;
      //std::cout << " linear break " << std::endl;
      break;
    }
    N = localspace.GetValOfShapeFunctionsAt(xietaphi).transpose().row(0);
    f = N*xcoor - targetPoint;
    //std::cout << "********************N*********************" << std::endl;
    //std::cout <<  N << std::endl;
    //std::cout << "********************N*xcoor *********************" << std::endl;
    //std::cout <<  N*xcoor  << std::endl;
    //std::cout << "********************targetPoint *********************" << std::endl;
    //std::cout <<  targetPoint  << std::endl;
    //std::cout << "********************xietaphi*********************" << std::endl;
    //std::cout <<  xietaphi<< std::endl;
    //std::cout << "********************dxietaphi*********************" << std::endl;
    //std::cout <<  dxietaphi<< std::endl;
    //std::cout << "***********************f******************" << std::endl;
    //std::cout <<  f << std::endl;
    //std::cout << "*****************************************" << std::endl;
    if( dxietaphi.squaredNorm() < 1e-3 && f.squaredNorm() < 1e-3 ){
      //std::cout << " tolerance break " << std::endl;
      break;
    }
  }
  if( x == 15){
    result.error = true;
    result.inside = false;
    //throw "error";
    return result;
  }
  result.error = false;
  //std::cout << " xietaphi "<< xietaphi << std::endl;
  result.distv = f;
  result.bary = xietaphi;
  result.baryClamped = ClampParamCoordinates(localspace.geoSupport, xietaphi);
  for(int p = elemDim; p< 3; ++p){
    result.baryClamped(p,0) = 0;
    result.bary(p,0) = 0;
  }
  // we need a tolerance here because in the case the source and the targe mesh are the same
  // we want to clasify all the point inside
  result.inside = (result.baryClamped-result.bary).squaredNorm() < 1e-5;

  return result;
}

};  // namespace BasicTools
