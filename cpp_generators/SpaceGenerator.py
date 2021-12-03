# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
#
from sympy import cse, ccode

import BasicTools.FE.Spaces.FESpaces as FES

from cpp_generators.Tools import PrintHeader, PrintToFile

def GetGeneratedFiles(prefix = "cpp_src"):
    return (prefix + "FE/GeneratedSpaces.cpp",)

def Generate(prefix = "cpp_src"):

    hfilename = prefix + "/FE/GeneratedSpaces.h"
    with open(hfilename,"w", encoding="utf8") as hfile:
        PrintHeader(hfile)
        PrintToFile(hfile,"#include <FE/Space.h>")
        PrintToFile(hfile,"const Space& GetSpaceFor(const std::string& spacename);")
        PrintToFile(hfile,"const std::vector<std::string> GetAvailableSpaces();")


    filename = prefix + "/FE/GeneratedSpaces.cpp"
    spaces = [("LagrangeSpaceGeo", FES.LagrangeSpaceGeo),
              ("ConstantSpaceGlobal", FES.ConstantSpaceGlobal),
              ("LagrangeSpaceP0", FES.LagrangeSpaceP0),
              ("LagrangeSpaceP1", FES.LagrangeSpaceP1),
              ("LagrangeSpaceP2", FES.LagrangeSpaceP2)]

    with open(filename,"w", encoding="utf8") as cppfile:
        PrintHeader(cppfile)
        PrintToFile(cppfile,"#include <memory>")
        PrintToFile(cppfile,"#include <stdexcept>")
        PrintToFile(cppfile,"#include <cmath>")
        PrintToFile(cppfile,"#include <LinAlg/EigenTypes.h>")
        PrintToFile(cppfile,"#include <FE/Space.h>")
        PrintToFile(cppfile,"using std::pow;")
        PrintToFile(cppfile,"using namespace BasicTools;")
        PrintToFile(cppfile,"")

        PrintToFile(cppfile,"namespace BasicTools{")

        for FEspn, FEspd in spaces:
            for spn in FEspd:
                spd = FEspd[spn]
                spd.Create()
                NOSF = spd.GetNumberOfShapeFunctions()
                nbDim = spd.GetDimensionality()

                cppfile.write(f"""
const int {FEspn}_{spn}_GetNumberOfShapeFunctions(){{ return {NOSF}; }};
const int {FEspn}_{spn}_GetDimensionality(){{ return {nbDim}; }};""")
                cppfile.write(f"""
MatrixDDD {FEspn}_{spn}_GetShapeFunctions(const double& phi, const double&  xi, const double&  eta  ){{
    MatrixDDD res = MatrixDDD::Zero({NOSF},1);
""")
                repl, redu = cse(spd.symN)
                for i, v in repl:
                    cppfile.write(f"    const double {i} = " +  ccode(v) + ";\n" )
                for ns in range(NOSF):
                    cppfile.write( f"    res.coeffRef({ns},1) =  "+ ccode(redu[0][ns]) )
                    cppfile.write(";\n")
                cppfile.write("""    return res;
};
""")
                cppfile.write(f"""
MatrixDDD {FEspn}_{spn}_GetShapeFunctionsDer(const double& phi, const double&  xi, const double&  eta  ){{
    MatrixDDD res = MatrixDDD::Zero({NOSF},3);
""")
                repl, redu = cse(spd.symdNdxi)
                for i, v in repl:
                    cppfile.write(f"    const double {i} = " +  ccode(v) + ";\n" )

                for j in range(NOSF):
                    for i in range(nbDim):
                        text = ccode(redu[0][i,j])
                        cppfile.write( f"    res.coeffRef({i},{j}) =  "+ text + ";\n")

                cppfile.write("""    return res;\n};
""")

        cppfile.write("""

std::map<std::string,Space> GetBasicSpaceAlmanac(){
std::map<std::string,Space> SpacesAlmanac;

\n""")
        for FEspn, FEspd in spaces:
            PrintToFile(cppfile,f"{{// working on space {FEspn}")
            PrintToFile(cppfile,"    Space localsp;")

            for spn in FEspd:
                spd = FEspd[spn]
                spd.Create()

                PrintToFile(cppfile,f"    {{// working on space: {spn}")
                PrintToFile(cppfile,"        ElementSpace fesp;")
                NOSF = spd.GetNumberOfShapeFunctions()

                for nsf in range(NOSF):
                    on,idxI,idxII = spd.dofAttachments[nsf]
                    if idxI is None :
                        idxI = -1
                    if idxII is None :
                        idxII = -1
                    if on == "F2" :
                        on = "E"
                    PrintToFile(cppfile,f"        fesp.AppendDofAttachement('{on[0]}', {idxI},{idxII});")
                PrintToFile(cppfile,f"        fesp.SFV = &{FEspn}_{spn}_GetShapeFunctions;")
                PrintToFile(cppfile,f"        fesp.SFDV = &{FEspn}_{spn}_GetShapeFunctionsDer;")
                PrintToFile(cppfile,f'        localsp.storage["{spn}"] = fesp;')
                PrintToFile(cppfile,f"    }};// end of space: {spn}")
            PrintToFile(cppfile,f'    SpacesAlmanac["{FEspn}"] = localsp;')
            PrintToFile(cppfile,f"}};// end of space: {FEspn}")
        PrintToFile(cppfile,"return SpacesAlmanac; \n}")


        PrintToFile(cppfile,"std::map<std::string,Space> SpacesAlmanacI = GetBasicSpaceAlmanac();")

        PrintToFile(cppfile,"""const Space& GetSpaceFor(const std::string& spacename){
    return SpacesAlmanacI[spacename];
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
