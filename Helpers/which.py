# -*- coding: utf-8 -*-
# From https://stackoverflow.com/questions/377017/test-if-executable-exists-in-python


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

def CheckIntegrity():
    return "OK"

if __name__ == '__main__':
    print((CheckIntegrity()))# pragma: no cover
