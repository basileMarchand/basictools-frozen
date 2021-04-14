//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//

#include <Containers/ElementFilter.h>

const MatrixID1 ElementFilterBase::GetIdsToTreatComplementaty( NativeUnstructuredMesh& mesh, const std::string& elemtype)  {
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
}
//
const MatrixID1 ElementFilterIntersection::GetIdsToTreat(NativeUnstructuredMesh& mesh, const std::string& elemtype)  {
    MatrixID1 res = storage[0]->GetIdsToTreat(mesh, elemtype);
    for (unsigned int i = 0; i< this->storage.size(); ++i){
        res = Intersect1D(res,storage[i]->GetIdsToTreat(mesh, elemtype));
    }
    return res;
};
//
