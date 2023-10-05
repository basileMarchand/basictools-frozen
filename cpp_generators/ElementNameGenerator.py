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

const GeoSupport GeoNA("NA",-1);
const GeoSupport GeoPoint("point",0);
const GeoSupport GeoBar("bar"  ,1);
const GeoSupport GeoTri("tri"  ,2);
const GeoSupport GeoQuad("quad" ,2);
const GeoSupport GeoTet("tet"  ,3);
const GeoSupport GeoPyr("pyr"  ,3);
const GeoSupport GeoWed("wed"  ,3);
const GeoSupport GeoHex("hex"  ,3 );

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
                #PrintToFile(cppFile,f"""        ElementNames["{elementType}"].faces.push_back(std::pair<ElementInfo,MatrixID1>(ElementNames["{e}"],faces) );""" )
                PrintToFile(cppFile,f"""        ElementNames["{elementType}"].faces.emplace_back(ElementNames["{e}"], faces);""" )
                PrintToFile(cppFile,"    }")
            for e,n in ei.faces2:
                PrintToFile(cppFile,"    {")
                PrintToFile(cppFile,"        MatrixID1 faces2;")
                PrintFillVMatrix(cppFile,"""        faces2""",n)
                PrintToFile(cppFile,f"""        ElementNames["{elementType}"].faces2.emplace_back(ElementNames["{e}"], faces2);""" )
                PrintToFile(cppFile,"    }")
            for e,n in ei.faces3:
                PrintToFile(cppFile,"    {")
                PrintToFile(cppFile,"        MatrixID1 faces3;")
                PrintFillVMatrix(cppFile,"""        faces3""",n)
                PrintToFile(cppFile,f"""        ElementNames["{elementType}"].faces3.emplace_back(ElementNames["{e}"], faces3);""" )
                PrintToFile(cppFile,"    }")


        PrintToFile(cppFile,"""
    return ElementNames;
};

std::map<std::string,ElementInfo> InitElementNames() ;
std::map<std::string,ElementInfo> ElementNames = InitElementNames();

}; // namespace BasicTools""")

if __name__ == '__main__':# pragma: no cover
    import sys
    sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
    if len(sys.argv) == 3  and sys.argv[1] == "-n":
        print(GetGeneratedFiles()[int(sys.argv[2])])
    elif len(sys.argv) == 3  and sys.argv[1] == "-g":
        Generate(sys.argv[2])
