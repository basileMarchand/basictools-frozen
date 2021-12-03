import numpy as np

cpp_header = """//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//

// this file was generated all changes will be lost!!!!!!!!!
"""

def PrintToFile(fp,message,rc = True):
    fp.write(message)
    if rc:
        fp.write("\n")

def PrintHeader(fp):
    PrintToFile(fp,cpp_header)

def PrintHeaderH(fp):
    PrintHeader(fp)
    PrintToFile(fp,"#pragma once")

def PrintFillVMatrix(fp, prefix, data):
    PrintToFile(fp,f"""{prefix}.resize({len(data)},1);""")
    PrintToFile(fp,f"""{prefix}  << """ + ", ".join(map(str,data)) + ";")

def PrintFillMatrix(fp, prefix, data):
    data = np.asarray(data)
    if len(data.shape) == 1:
        PrintToFile(fp,f"""{prefix}.resize(1,{len(data)});""")
        PrintToFile(fp,f"""{prefix}  << """ + ", ".join(map(str,data)) + ";")
    else:
        PrintToFile(fp,f"""{prefix}.resize{str(data.shape)};""")
        for i in range(data.shape[0]):
            PrintToFile(fp,f"""{prefix}.row({i})  << """ + ", ".join(map(str,data[i,:])) + ";")

def PrintBool(fp, prefix, data):
    strdata = 'true' if data else 'false'
    PrintToFile(fp,prefix + f"= {strdata};")