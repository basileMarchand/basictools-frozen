//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//
#include <Containers/UnstructuredMeshTools.h>
//#include <LinAlg/EigenTypes.h>
#include <LinAlg/EigenTools.h>

namespace BasicTools
{
    
MatrixDDD GetElementsCenters(const MapMatrixDDD& nodes, const ElementsContainer& elements){

     auto& connectivity = elements.GetConnectivityMatrix();
     MatrixDDD res(elements.GetNumberOfElements(), nodes.cols());
     res.fill(0.0);
     for(Eigen::Index i=0; i <nodes.cols(); ++i){
         MatrixDDD a = indexingi(nodes,connectivity,i);

         MatrixDD1 b= a.rowwise().sum().eval().col(0);

         //MatrixDD1 c = (b.col(0)+ res.col(i)).eval();
         res.col(i) = (b.col(0)+ res.col(i)).eval();
         //res.col(i) += indexingi(nodes,connectivity,i).rowwise().sum();
     }
     res /= connectivity.cols();
     return res;

};

} // namespace BasicTools
