# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
from typing import TextIO
from BasicTools.NumpyDefs import ArrayLike
cpp_header = """//
// This file is subject to the terms and conditions defined in
// file 'LICENSE.txt', which is part of this source code package.
//

// this file was generated all changes will be lost!!!!!!!!!
"""

def PrintToFile(fp:TextIO, message:str, rc:bool = True):
    """Write text to a file

    Parameters
    ----------
    fp : TextIO
        the file descriptor
    message : str
        the message to write
    rc : bool, optional
        if a "\n" must be added at the end, by default True
    """
    fp.write(message)
    if rc:
        fp.write("\n")

def PrintHeader(fp:TextIO):
    """Print the BasicTools standard .cpp header

    Parameters
    ----------
    fp : TextIO
        the file descriptor
    """
    PrintToFile(fp,cpp_header)

def PrintHeaderH(fp:TextIO):
    """Print the BasicTools standard .h header.
    This include the "#pragma once" line

    Parameters
    ----------
    fp : TextIO
        the file descriptor
    """
    PrintHeader(fp)
    PrintToFile(fp,"#pragma once")

def PrintFillVMatrix(fp:TextIO, prefix:str, data:ArrayLike):
    """ Write vector in vertical

    Parameters
    ----------
    fp : TextIO
        the file descriptor
    prefix : str
        the name of the variable (this can include indentation spaces)
    data : ArrayLike
        the value of the variable
    """
    PrintToFile(fp,f"""{prefix}.resize({len(data)},1);""")
    PrintToFile(fp,f"""{prefix}  << """ + ", ".join(map(str,data)) + ";")

def PrintFillMatrix(fp:TextIO, prefix:str, data:ArrayLike):
    """ Write vector in horizontal format

    Parameters
    ----------
    fp : TextIO
        the file descriptor
    prefix : str
        the name of the variable (this can include indentation spaces)
    data : ArrayLike
        the value of the variable
    """
    import numpy as np
    data = np.asarray(data)
    if len(data.shape) == 1:
        PrintToFile(fp,f"""{prefix}.resize(1,{len(data)});""")
        PrintToFile(fp,f"""{prefix}  << """ + ", ".join(map(str,data)) + ";")
    else:
        PrintToFile(fp,f"""{prefix}.resize{str(data.shape)};""")
        for i in range(data.shape[0]):
            PrintToFile(fp,f"""{prefix}.row({i})  << """ + ", ".join(map(str,data[i,:])) + ";")

def PrintBool(fp:TextIO, prefix:str, data:bool):
    """Helper function to set a variable to a boolean value

    Parameters
    ----------
    fp : TextIO
        the file descriptor
    prefix : str
        the name of the variable (this can include indentation spaces)
    data : ArrayLike
        the value of the variable
    """
    strData = 'true' if data else 'false'
    PrintToFile(fp,f"{prefix}= {strData};")