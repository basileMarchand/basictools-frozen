# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
#
import os
from typing import Tuple

def GetGeneratedFiles(prefix: str = "cpp_src") -> Tuple[str]:
    """ Get the list of generated files for this generator

    Parameters
    ----------
    prefix : str, optional
        prefix for the generated files, by default "cpp_src/"

    Returns
    -------
    Tuple[str]
        list of generated files for this generator
    """
    return ( os.path.join(prefix ,"Containers", "GeneratedElementNames.cpp"), )

def Generate(prefix:str = "cpp_src") :
    """ Run the generation of cpp file using the prefix.
    the file created: prefix + "Containers/GeneratedElementNames.cpp"
    This file contains all the element related data: name, dimensionality,
    faces...

    Parameters
    ----------
    prefix : str, optional
        prefix for the generated files, by default "cpp_src"


    """
    from BasicTools.Containers.ElementNames import ElementsInfo
    from cpp_generators.Tools import PrintHeader, PrintToFile, PrintFillMatrix, PrintBool, PrintFillVMatrix
    filename = GetGeneratedFiles(prefix)[0]

    with open(filename,"w", encoding="utf8") as cppFile:
        PrintHeader(cppFile)
        cppFile.write("""
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

std::map<std::string,ElementInfo> InitElementNames() {
    std::map<std::string,ElementInfo> ElementNames;
""")
        for elementType,ei in ElementsInfo.items() :

            PrintToFile(cppFile,f"""
    ElementNames["{elementType}"] = ElementInfo();
    ElementNames["{elementType}"].numberOfNodes = {ei.numberOfNodes};
    ElementNames["{elementType}"].geoSupport = Geo{ei.geoSupport.name.capitalize()};
    ElementNames["{elementType}"].name = "{elementType}";""")

            PrintFillMatrix(cppFile,f"""    ElementNames["{elementType}"].mirrorPermutation""",ei.mirrorPermutation)
            PrintBool(cppFile,f"""    ElementNames["{elementType}"].linear """, ei.linear)
            PrintToFile(cppFile,f"""    ElementNames["{elementType}"].degree = {ei.degree}; """ )
            for e,n in ei.faces:
                PrintToFile(cppFile,"    {")
                PrintToFile(cppFile,"        MatrixID1 faces;")
                PrintFillVMatrix(cppFile,"""        faces""",n)
                PrintToFile(cppFile,f"""        ElementNames["{elementType}"].faces.push_back(std::pair<ElementInfo,MatrixID1>(ElementNames["{e}"],faces) );""" )
                PrintToFile(cppFile,"    }")
            for e,n in ei.faces2:
                PrintToFile(cppFile,"    {")
                PrintToFile(cppFile,"        MatrixID1 faces2;")
                PrintFillVMatrix(cppFile,"""        faces2""",n)
                PrintToFile(cppFile,f"""        ElementNames["{elementType}"].faces2.push_back(std::pair<ElementInfo,MatrixID1>(ElementNames["{e}"],faces2) );""" )
                PrintToFile(cppFile,"    }")
        PrintToFile(cppFile,"""
    return ElementNames;
};

std::map<std::string,ElementInfo> InitElementNames() ;
std::map<std::string,ElementInfo> ElementNames = InitElementNames();

}; // namespace BasicTools""")
