# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
import os
from typing import Tuple

def GetGeneratedFiles(prefix:str  = "cpp_src") -> Tuple[str]:
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
    return ( os.path.join(prefix ,"FE", "GeneratedIntegrationsRules.cpp"), )

def Generate(prefix:str = "cpp_src"):
    """ Run the generation of cpp file using the prefix.
    the file created: prefix + "Containers/GeneratedElementNames.cpp"
    This file contains all the integration rules defined in the python
    files.

    Parameters
    ----------
    prefix : str, optional
        prefix for the generated files, by default "cpp_src"


    """
    from cpp_generators.Tools import PrintHeader, PrintToFile, PrintFillMatrix
    import BasicTools.FE.IntegrationsRules as IR

    """Run the generation of cpp file using the prefix"""
    filename = GetGeneratedFiles(prefix)[0]

    with open(filename, 'w', encoding="utf8") as cppFile:
        PrintHeader(cppFile)
        PrintToFile(cppFile,"#include <map>")
        PrintToFile(cppFile,"#include <stdexcept>")
        PrintToFile(cppFile,"#include <cmath>")
        PrintToFile(cppFile,"#include <LinAlg/EigenTypes.h>")
        PrintToFile(cppFile,"#include <FE/Space.h>")
        PrintToFile(cppFile,"#include <FE/IntegrationRule.h>")
        PrintToFile(cppFile,"using std::pow;")
        PrintToFile(cppFile,"using namespace BasicTools;")
        PrintToFile(cppFile,"""
namespace BasicTools {

std::map<std::string,SpaceIntegrationRule>  GetPythonDefinedIntegrationRules(){
    std::map<std::string,SpaceIntegrationRule> res;
""")
        for k,v in IR.IntegrationRulesAlmanac.items():

            PrintToFile(cppFile,f"""    {{ // working on {k}
            SpaceIntegrationRule ir; """)
            for k2,v2 in v.items():
                PrintToFile(cppFile,f"""        {{ // working on {k2}
                IntegrationRule eir;""")
                PrintFillMatrix(cppFile,12*" "+"eir.p", v2[0])
                PrintFillMatrix(cppFile,12*" "+"eir.w", v2[1])
                PrintToFile(cppFile,12*" "+f"""ir.storage["{k2}"] = eir;""" )
                PrintToFile(cppFile,8*" "+"}")
            PrintToFile(cppFile,8*" "+f"""res["{k}"] = ir;""" )
            PrintToFile(cppFile,8*" "+"}")
        PrintToFile(cppFile,"""
    return res;
};

std::map<std::string,SpaceIntegrationRule> IntegrationRulesAlmanac = GetPythonDefinedIntegrationRules();
};// BasicTools namespace
""")

if __name__ == '__main__':# pragma: no cover
    import sys
    sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
    if len(sys.argv) == 3  and sys.argv[1] == "-n":
        print(GetGeneratedFiles()[int(sys.argv[2])])
    elif len(sys.argv) == 3  and sys.argv[1] == "-g":
        Generate(sys.argv[2])

