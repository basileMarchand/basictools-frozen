//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//
#pragma once

#include <utility>
#include <string>
#include <LinAlg/EigenTypes.h>

class GeoSupport {
public:
    std::string name;
    int dimensionality;
    GeoSupport(const std::string& name, const int& dimensionality): name(name), dimensionality(dimensionality) {}
    std::string ToStr() const {
        return std::string("GeoSuport( " + this->name + ")");
    }
};

GeoSupport GeoNA("NA",-1);
GeoSupport GeoPoint("point",0);
GeoSupport GeoBar("bar"  ,1);
GeoSupport GeoTri("tri"  ,2);
GeoSupport GeoQuad("quad" ,2);
GeoSupport GeoTet("tet"  ,3);
GeoSupport GeoPyr("pyr"  ,3);
GeoSupport GeoWed("wed"  ,3);
GeoSupport GeoHex("hex"  ,3 );

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

const std::string Point_1 = "point1";
const std::string Bar_2 = "bar2";
const std::string Bar_3 = "bar3";
const std::string Triangle_6 = "tri6";
const std::string Quadrangle_8 = "quad8";
const std::string Quadrangle_9 = "quad9";
const std::string Tetrahedron_4 = "tet4";
const std::string Pyramid_5 = "pyr5";
const std::string Wedge_6 = "wed6";
const std::string Hexaedron_8 = "hex8";
const std::string Tetrahedron_10 = "tet10";
const std::string Pyramid_13 = "pyr13";
const std::string Wedge_15 = "wed15";
const std::string Wedge_18 = "wed18";
const std::string Hexaedron_20 = "hex20";
const std::string Hexaedron_27 = "hex27";

std::map<std::string,ElementInfo> InitElementNames() ;

std::map<std::string,ElementInfo> ElementNames = InitElementNames();
