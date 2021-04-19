//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//
#pragma once

#include <utility>
#include <string>
#include <LinAlg/EigenTypes.h>

namespace BasicTools
{
    
class GeoSupport {
public:
    std::string name;
    int dimensionality;
    GeoSupport(const std::string& name, const int& dimensionality): name(name), dimensionality(dimensionality) {}
    std::string ToStr() const {
        return std::string("GeoSuport( " + this->name + ")");
    }
};

extern GeoSupport GeoNA;
extern GeoSupport GeoPoint;
extern GeoSupport GeoBar;
extern GeoSupport GeoTri;
extern GeoSupport GeoQuad;
extern GeoSupport GeoTet;
extern GeoSupport GeoPyr;
extern GeoSupport GeoWed;
extern GeoSupport GeoHex;

class ElementInfo{
public:
    int numberOfNodes;
    GeoSupport geoSupport;
    std::string name;
    MatrixID1 mirrorPermutation;
    bool linear;
    std::vector<std::pair<ElementInfo,MatrixID1> > faces;
    std::vector<std::pair<ElementInfo,MatrixID1> > faces2;
    ElementInfo(): geoSupport("NA",-1){}
    int dimension() const { return this->geoSupport.dimensionality; }
};

extern const std::string Point_1;
extern const std::string Bar_2;
extern const std::string Bar_3;
extern const std::string Triangle_6;
extern const std::string Quadrangle_8;
extern const std::string Quadrangle_9;
extern const std::string Tetrahedron_4;
extern const std::string Pyramid_5;
extern const std::string Wedge_6;
extern const std::string Hexaedron_8;
extern const std::string Tetrahedron_10;
extern const std::string Pyramid_13;
extern const std::string Wedge_15;
extern const std::string Wedge_18;
extern const std::string Hexaedron_20;
extern const std::string Hexaedron_27;

extern std::map<std::string,ElementInfo> ElementNames;

} // namespace BasicTools
