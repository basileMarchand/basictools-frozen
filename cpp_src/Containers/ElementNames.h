//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//
#pragma once

#include <cassert>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include <LinAlg/EigenTypes.h>

namespace BasicTools {

struct GeoSupport {
    std::string name{"NA"};
    int dimensionality{-1};

    GeoSupport() = default;
    GeoSupport(const std::string& name, const int& dimensionality) : name(name), dimensionality(dimensionality) {}

    std::string ToStr() const {
        return std::string("GeoSuport( " + this->name + ")");
    }
};

extern const GeoSupport GeoNA;
extern const GeoSupport GeoPoint;
extern const GeoSupport GeoBar;
extern const GeoSupport GeoTri;
extern const GeoSupport GeoQuad;
extern const GeoSupport GeoTet;
extern const GeoSupport GeoPyr;
extern const GeoSupport GeoWed;
extern const GeoSupport GeoHex;

struct ElementInfo {
    int numberOfNodes;
    GeoSupport geoSupport;
    std::string name;
    MatrixI1D mirrorPermutation;
    bool linear;
    int degree;
    std::vector<std::pair<ElementInfo,MatrixID1>> faces;
    std::vector<std::pair<ElementInfo,MatrixID1>> faces2;
    std::vector<std::pair<ElementInfo,MatrixID1>> faces3;

    ElementInfo() = default;

    int dimension() const { return this->geoSupport.dimensionality; }
    const std::vector<std::pair<ElementInfo, MatrixID1>>& GetFacesLevel(int level) {
        assert(level > 0);
        assert(level < 3);
        if (level==1) return faces;
        if (level==2) return faces2;
        if (level==3) return faces3;
        throw;
    }
};

inline const std::string Point_1 = "point1";
inline const std::string Bar_2 = "bar2";
inline const std::string Bar_3 = "bar3";
inline const std::string Triangle_3 = "tri3";
inline const std::string Triangle_6 = "tri6";
inline const std::string Quadrangle_4 = "quad4";
inline const std::string Quadrangle_8 = "quad8";
inline const std::string Quadrangle_9 = "quad9";
inline const std::string Tetrahedron_4 = "tet4";
inline const std::string Pyramid_5 = "pyr5";
inline const std::string Wedge_6 = "wed6";
inline const std::string Hexaedron_8 = "hex8";
inline const std::string Tetrahedron_10 = "tet10";
inline const std::string Pyramid_13 = "pyr13";
inline const std::string Wedge_15 = "wed15";
inline const std::string Wedge_18 = "wed18";
inline const std::string Hexaedron_20 = "hex20";
inline const std::string Hexaedron_27 = "hex27";

extern std::map<std::string,ElementInfo> ElementNames;

} // namespace BasicTools
