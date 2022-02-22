# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

"""Class to check if a executable is in the current path

insperated from:
    https://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
"""
def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
            if os.name == "nt":# Windows
                try:
                    from win32api import FindExecutable, GetLongPathName
                    _, executable = FindExecutable(program)
                    if os.path.isfile(executable):
                        return executable
                except:
                    pass

    return None

def CheckIntegrity():
    print(which("ls"))
    print(which("dir"))
    return "OK"

if __name__ == '__main__':
    print((CheckIntegrity()))# pragma: no cover
