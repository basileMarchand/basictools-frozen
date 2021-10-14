import BasicTools.FE.IntegrationsRules as IR

from cpp_generators.aux import PrintHeader, PrintToFile, PrintFillMatrix

def GetGeneratedFiles(prefix = "cpp_src/"):
    return (prefix + "FE/GeneratedIntegrationsRules.cpp",)

def Generate(prefix = "cpp_src/"):
    filename = prefix +  "FE/GeneratedIntegrationsRules.cpp"


    cppfile = open(filename,"w")
    PrintHeader(cppfile)
    PrintToFile(cppfile,"#include <map>")
    PrintToFile(cppfile,"#include <stdexcept>")
    PrintToFile(cppfile,"#include <cmath>")
    PrintToFile(cppfile,"#include <LinAlg/EigenTypes.h>")
    PrintToFile(cppfile,"#include <FE/Space.h>")
    PrintToFile(cppfile,"#include <FE/IntegrationRule.h>")
    PrintToFile(cppfile,"using std::pow;")
    PrintToFile(cppfile,"using namespace BasicTools;")
    PrintToFile(cppfile,"""
namespace BasicTools { 

std::map<std::string,SpaceIntegrationRule>  GetPythonDefinedIntegrationRules(){
    std::map<std::string,SpaceIntegrationRule> res;
""")
    for k,v in IR.IntegrationRulesAlmanac.items():

        PrintToFile(cppfile,f"""    {{ // working on {k} 
        SpaceIntegrationRule ir; """)
        for k2,v2 in v.items():
            PrintToFile(cppfile,f"""        {{ // working on {k2} 
            IntegrationRule eir;""")
            PrintFillMatrix(cppfile,12*" "+"eir.p", v2[0])
            PrintFillMatrix(cppfile,12*" "+"eir.w", v2[1])
            PrintToFile(cppfile,12*" "+f"""ir.storage["{k2}"] = eir;""" )
            PrintToFile(cppfile,8*" "+"}")
        PrintToFile(cppfile,8*" "+f"""res["{k}"] = ir;""" )
        PrintToFile(cppfile,8*" "+"}")
    PrintToFile(cppfile,"""
    return res;
};

std::map<std::string,SpaceIntegrationRule> IntegrationRulesAlmanac = GetPythonDefinedIntegrationRules();
};// BasicTools namespace
""")

if __name__ == '__main__':# pragma: no cover
    Generate()

