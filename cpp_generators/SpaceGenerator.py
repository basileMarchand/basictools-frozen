# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
#

def GetGeneratedFiles(prefix = "cpp_src"):
    return (prefix + "FE/GeneratedSpaces.cpp",)

def Generate(prefix = "cpp_src"):
    from cpp_generators.Tools import PrintHeader, PrintToFile
    from sympy import cse, ccode
    import BasicTools.FE.Spaces.FESpaces as FES

    hFileName = prefix + "/FE/GeneratedSpaces.h"
    with open(hFileName,"w", encoding="utf8") as hfile:
        PrintHeader(hfile)
        PrintToFile(hfile,"#include <FE/Space.h>")
        PrintToFile(hfile,"const Space& GetSpaceFor(const std::string& spaceName);")
        PrintToFile(hfile,"const std::vector<std::string> GetAvailableSpaces();")


    filename = prefix + "/FE/GeneratedSpaces.cpp"
    spaces = [("LagrangeSpaceGeo", FES.LagrangeSpaceGeo),
              ("ConstantSpaceGlobal", FES.ConstantSpaceGlobal),
              ("LagrangeSpaceP0", FES.LagrangeSpaceP0),
              ("LagrangeSpaceP1", FES.LagrangeSpaceP1),
              ("LagrangeSpaceP2", FES.LagrangeSpaceP2)]

    with open(filename,"w", encoding="utf8") as cppFile:
        PrintHeader(cppFile)
        PrintToFile(cppFile,"#include <memory>")
        PrintToFile(cppFile,"#include <stdexcept>")
        PrintToFile(cppFile,"#include <cmath>")
        PrintToFile(cppFile,"#include <LinAlg/EigenTypes.h>")
        PrintToFile(cppFile,"#include <FE/Space.h>")
        PrintToFile(cppFile,"using std::pow;")
        PrintToFile(cppFile,"using namespace BasicTools;")
        PrintToFile(cppFile,"")

        PrintToFile(cppFile,"namespace BasicTools{")

        for FESpaceName, FEspace in spaces:
            for spn in FEspace:
                spd = FEspace[spn]
                spd.Create()
                numberOfShapeFunctions = spd.GetNumberOfShapeFunctions()
                nbDim = spd.GetDimensionality()

                cppFile.write(f"""
const int {FESpaceName}_{spn}_GetNumberOfShapeFunctions(){{ return {numberOfShapeFunctions}; }};
const int {FESpaceName}_{spn}_GetDimensionality(){{ return {nbDim}; }};""")
                cppFile.write(f"""
MatrixDDD {FESpaceName}_{spn}_GetShapeFunctions(const double& phi, const double&  xi, const double&  eta  ){{
    MatrixDDD res = MatrixDDD::Zero({numberOfShapeFunctions},1);
""")
                replacements, reducedExprs  = cse(spd.symN)
                for i, v in replacements:
                    cppFile.write(f"    const double {i} = " +  ccode(v) + ";\n" )
                for ns in range(numberOfShapeFunctions):
                    cppFile.write( f"    res.coeffRef({ns},1) =  "+ ccode(reducedExprs [0][ns]) )
                    cppFile.write(";\n")
                cppFile.write("""    return res;
};
""")
                cppFile.write(f"""
MatrixDDD {FESpaceName}_{spn}_GetShapeFunctionsDer(const double& phi, const double&  xi, const double&  eta  ){{
    MatrixDDD res = MatrixDDD::Zero({numberOfShapeFunctions},3);
""")
                replacements, reducedExprs  = cse(spd.symdNdxi)
                for i, v in replacements:
                    cppFile.write(f"    const double {i} = " +  ccode(v) + ";\n" )

                for j in range(numberOfShapeFunctions):
                    for i in range(nbDim):
                        text = ccode(reducedExprs [0][i,j])
                        cppFile.write( f"    res.coeffRef({i},{j}) =  "+ text + ";\n")

                cppFile.write("""    return res;\n};
""")

        cppFile.write("""

std::map<std::string,Space> GetBasicSpaceAlmanac(){
std::map<std::string,Space> SpacesAlmanac;

\n""")
        for FESpaceName, FEspace in spaces:
            PrintToFile(cppFile,f"{{// working on space {FESpaceName}")
            PrintToFile(cppFile,"    Space localsp;")

            for spn in FEspace:
                spd = FEspace[spn]
                spd.Create()

                PrintToFile(cppFile,f"    {{// working on space: {spn}")
                PrintToFile(cppFile,"        ElementSpace fesp;")
                numberOfShapeFunctions = spd.GetNumberOfShapeFunctions()

                for nsf in range(numberOfShapeFunctions):
                    on,idxI,idxII = spd.dofAttachments[nsf]
                    if idxI is None :
                        idxI = -1
                    if idxII is None :
                        idxII = -1
                    if on == "F2" :
                        on = "E"
                    PrintToFile(cppFile,f"        fesp.AppendDofAttachment('{on[0]}', {idxI},{idxII});")
                PrintToFile(cppFile,f"        fesp.SFV = &{FESpaceName}_{spn}_GetShapeFunctions;")
                PrintToFile(cppFile,f"        fesp.SFDV = &{FESpaceName}_{spn}_GetShapeFunctionsDer;")
                PrintToFile(cppFile,f'        localsp.storage["{spn}"] = fesp;')
                PrintToFile(cppFile,f"    }};// end of space: {spn}")
            PrintToFile(cppFile,f'    SpacesAlmanac["{FESpaceName}"] = localsp;')
            PrintToFile(cppFile,f"}};// end of space: {FESpaceName}")
        PrintToFile(cppFile,"return SpacesAlmanac; \n}")


        PrintToFile(cppFile,"std::map<std::string,Space> SpacesAlmanacI = GetBasicSpaceAlmanac();")

        PrintToFile(cppFile,"""const Space& GetSpaceFor(const std::string& spaceName){
    return SpacesAlmanacI[spaceName];
};

const std::vector<std::string> GetAvailableSpaces(){
    std::vector<std::string> res ;
    for(auto const& item: SpacesAlmanacI ) res.push_back(item.first);
    return res;
}
}; // BasicTools Namespace
""")

def CheckIntegrity(GUI=False):
    from BasicTools.Helpers.Tests import TestTempDir
    Generate(prefix = TestTempDir.GetTempPath())
    return "ok"

if __name__ == '__main__':# pragma: no cover
    print(CheckIntegrity(GUI=True))
