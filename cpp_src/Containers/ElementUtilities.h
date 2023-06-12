#pragma once

#include <LinAlg/EigenTypes.h>

namespace BasicTools {

MatrixDDD ClampParamCoordinates(const GeoSupport& geoSupport, MatrixD1D coord){

std::cout << "*-*-*-*------>"<< GeoBar.name  << "*"<< GeoBar.dimensionality  <<"<----------------------------" << coord<< std::endl;
std::cout << "coods original ------------------------>"<< geoSupport.name  << "*"<< geoSupport.dimensionality  <<"<----------------------------" << coord<< std::endl;
if (geoSupport.name == "point" ) {
    assert(false);
    return coord;
} else if (geoSupport.name == "bar" ) {
    coord(0,0) = coord(0,0)>1. ? 1. :  coord(0,0);
    coord(0,0) = coord(0,0)<0. ? 0. :  coord(0,0);
    std::cout << "coods clamped "<< coord << std::endl;
    return coord;
} else if (geoSupport.name == "tri" ) {
    assert(false);
    return coord;
} else if (geoSupport.name == "quad" ) {
    for(int i =0; i < 2;++i){
        coord(0,i) = coord(0,i)>1. ? 1. :  coord(0,i);
        coord(0,i) = coord(0,i)<0. ? 0. :  coord(0,i);
    }
    return coord;
} else if (geoSupport.name == "tet" ) {
    assert(false);
    return coord;
} else if (geoSupport.name == "pyr" ) {
    assert(false);
    return coord;
} else if (geoSupport.name == "wed" ) {
    assert(false);
    return coord;
} else if (geoSupport.name == "hex" ) {
    for(int i =0; i < 3;++i){
        coord(0,i) = coord(0,i)>1. ? 1. :  coord(0,i);
        coord(0,i) = coord(0,i)<0. ? 0. :  coord(0,i);
    }
    return coord;
} else if (geoSupport.name == "NA" ) {
    assert(false);
    return coord;
}
return coord;
}


};  // namespace BasicTools