#pragma once

#include <LinAlg/EigenTypes.h>

namespace BasicTools {


template<int dim>
void ClampComplex( MatrixD31& coord){
    for(int i =0; i < dim;++i){
        coord(i,0) = coord(i,0)>1. ? 1. :  coord(i,0);
        coord(i,0) = coord(i,0)<0. ? 0. :  coord(i,0);
    }
}

MatrixD31 ClampParamCoordinates(const GeoSupport& geoSupport, MatrixD31 coord){

//std::cout << "*-*-*-*------>"<< GeoBar.name  << "*"<< GeoBar.dimensionality  <<"<----------------------------==" << coord<< std::endl;
//std::cout << "coods original ------------------------>"<< geoSupport.name  << "*"<< geoSupport.dimensionality  <<"<----------------------------===" << coord<< std::endl;
if (geoSupport.name == "point" ) {
    return coord;
} else if (geoSupport.name == "bar" ) {
    ClampComplex<1>(coord);
    return coord;
} else if (geoSupport.name == "tri" ) {
    if(coord(0,0) + coord(1,0)  > 1 ){
        const double dif = coord(0,0)-coord(1,0);
        coord(0,0) = (1+dif)/2.;
        coord(1,0) = 1.-coord(0,0);
    }
    ClampComplex<2>(coord);
    return coord;
} else if (geoSupport.name == "quad" ) {
    ClampComplex<2>(coord);
    return coord;
} else if (geoSupport.name == "tet" ) {
    const double s = coord(0,0) + coord(1,0) + coord(2,0);
    if(s > 1){
       coord.array() += (1.-s)/3.;
    }
    ClampComplex<3>(coord);
    return coord;
} else if (geoSupport.name == "pyr" ) {
    throw "pyr not implemente in the field transfer";
    return coord;
} else if (geoSupport.name == "wed" ) {
     if(coord(0,0) + coord(1,0)  > 1 ){
        const double dif = coord(0,0)-coord(1,0);
        coord(0,0) = (1+dif)/2.;
        coord(1,0) = 1.-coord(0,0);
    }
    ClampComplex<3>(coord);
    return coord;
} else if (geoSupport.name == "hex" ) {
    ClampComplex<3>(coord);
    return coord;
} else if (geoSupport.name == "NA" ) {
    throw "NON NON NON!!!!";
    return coord;
}
throw "NON NON NON!!!!";
return coord;
}


};  // namespace BasicTools