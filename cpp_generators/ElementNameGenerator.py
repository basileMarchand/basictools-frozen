# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
#
from BasicTools.Containers.ElementNames import ElementsInfo
from cpp_generators.Tools import PrintHeader, PrintToFile, PrintFillMatrix, PrintBool, PrintFillVMatrix

def GetGeneratedFiles(prefix = "cpp_src/"):
    """Get the list of generated files for this generator"""
    return ( prefix+ "Containers/GeneratedElementNames.cpp",)

def Generate(prefix = "cpp_src/"):
    """Run the generation of cpp file using the prefix"""
    filename = prefix + "Containers/GeneratedElementNames.cpp"

    with open(filename,"w", encoding="utf8") as cppfile:
        PrintHeader(cppfile)
        cppfile.write("""
#include <Containers/ElementNames.h>
#include <map>

namespace BasicTools{

GeoSupport GeoNA("NA",-1);
GeoSupport GeoPoint("point",0);
GeoSupport GeoBar("bar"  ,1);
GeoSupport GeoTri("tri"  ,2);
GeoSupport GeoQuad("quad" ,2);
GeoSupport GeoTet("tet"  ,3);
GeoSupport GeoPyr("pyr"  ,3);
GeoSupport GeoWed("wed"  ,3);
GeoSupport GeoHex("hex"  ,3 );

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

std::map<std::string,ElementInfo> InitElementNames() {
    std::map<std::string,ElementInfo> ElementNames;
""")
        for elemtype,ei in ElementsInfo.items() :

            PrintToFile(cppfile,"""
    ElementNames["{elemtype}"] = ElementInfo();
    ElementNames["{elemtype}"].numberOfNodes = {non};
    ElementNames["{elemtype}"].geoSupport = Geo{geoname};
    ElementNames["{elemtype}"].name = "{elemtype}";""".format(elemtype=elemtype,non = ei.numberOfNodes,geoname=ei.geoSupport.name.capitalize()))

            PrintFillMatrix(cppfile,f"""    ElementNames["{elemtype}"].mirrorPermutation""",ei.mirrorPermutation)
            PrintBool(cppfile,f"""    ElementNames["{elemtype}"].linear """, ei.linear)
            PrintToFile(cppfile,f"""    ElementNames["{elemtype}"].degree = {ei.degree}; """ )
            for e,n in ei.faces:
                PrintToFile(cppfile,"    {")
                PrintToFile(cppfile,"        MatrixID1 faces;")
                PrintFillVMatrix(cppfile,"""        faces""",n)
                PrintToFile(cppfile,f"""        ElementNames["{elemtype}"].faces.push_back(std::pair<ElementInfo,MatrixID1>(ElementNames["{e}"],faces) );""" )
                PrintToFile(cppfile,"    }")
            for e,n in ei.faces2:
                PrintToFile(cppfile,"    {")
                PrintToFile(cppfile,"        MatrixID1 faces2;")
                PrintFillVMatrix(cppfile,"""        faces2""",n)
                PrintToFile(cppfile,f"""        ElementNames["{elemtype}"].faces2.push_back(std::pair<ElementInfo,MatrixID1>(ElementNames["{e}"],faces2) );""" )
                PrintToFile(cppfile,"    }")
        PrintToFile(cppfile,"""
    return ElementNames;
};

std::map<std::string,ElementInfo> InitElementNames() ;
std::map<std::string,ElementInfo> ElementNames = InitElementNames();

}; // namespace BasicTools""")
