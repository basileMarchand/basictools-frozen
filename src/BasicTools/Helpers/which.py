# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

from typing import Union

def which(program: str) -> Union[str,None]:
    """Check if an executable is reachable
    and return the real "executable" (with the .exe or .bat )

    inspiration from:
        https://stackoverflow.com/questions/377017/test-if-executable-exists-in-python

    Parameters
    ----------
    program : str
        the name of the programe to find, in windows the user can omit the .exe or the .bat
        extentions

    Returns
    -------
    Union[str,None]
        return the string of the real executable or None if the executable is not found
    """

    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)

    if fpath:
        # we have a path so we search for the exact filename
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
            if os.name == "nt":# Windows # pragma: no cover
                try:
                    from win32api import FindExecutable, GetLongPathName
                    _, executable = FindExecutable(program)
                except:
                    pass
                else:
                    if os.path.isfile(executable):
                        return executable

    return None

def CheckIntegrity():
    import os

    print(which("ls"))
    print(which("dir"))
    print(which("lls"))
    if os.name == "nt": # pragma: no cover
        pass
    else:
        print(which("/usr/bin/ls"))

    return "OK"

if __name__ == '__main__': # pragma: no cover
    print((CheckIntegrity()))
