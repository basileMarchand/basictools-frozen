# -*- coding: utf-8 -*-
"""Class to check if a executable is in the current path

insperated from:
    https://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
"""
__author__ = "Felipe Bordeu"

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

    return None