# -*- coding: utf-8 -*-
from __future__ import print_function

import sys
import os
"""
Send all the print statements to a file or to the sink os ...
"""


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

class PrintBypass():

    def __init__(self):
        self.stdout_ = sys.stdout #Keep track of the previous value.
        self.stderr_ = sys.stderr #Keep track of the previous value.

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.Restore()

    def ToSink(self):
        sys.stdout = open(os.devnull,"w")
        sys.stderr = open(os.devnull,"w")

    def ToDisk(self,filename, filenamecerr=None):
        sys.stdout = open(filename, 'w') # Something here that provides a write method.
        if filenamecerr is None:
            sys.stderr = sys.stdout
        else:
            sys.stderr = open(filenamecerr, 'w')

    def CopyOutputToDisk(self,filename, filenamecerr=None):
        class COTD():
            def __init__(self,buffertocopy,filename):
                self.Fdout = open(filename, 'w')
                self.sysout = buffertocopy

            def close(self):
                self.Fdout.close()
            def flush(self):
                self.Fdout.flush()
                self.sysout.flush()
            def write(self,data):
                self.Fdout.write(data)
                self.sysout.write(data)

        sys.stdout = COTD(self.stdout_, filename)
        if filenamecerr is None:
            sys.stderr = sys.stdout
        else:
            sys.stderr = COTD(self.stderr_,filenamecerr)

    def ToRedirect(self,cout_obj,cerr_obj=None):
        # obj must implement the functions: close(), flush(), write(data)
        sys.stdout = cout_obj
        if cerr_obj is not None:
            sys.stderr = cerr_obj

    def Restore(self):
        if sys.stdout is not self.stdout_:
            self.Print("Restore stdout")
            sys.stdout.flush()
            sys.stdout.close()
            sys.stdout = self.stdout_  # restore the previous stdout.
        if sys.stderr is not self.stderr_ :
            self.Print("Restore stdcerr")
            sys.stderr.flush()
            sys.stderr.close()
            sys.stderr = self.stderr_  # restore the previous stdout.

    def Print(self,text):
        """To print to the original cout"""
        self.stdout_.write(text+"\n")

    def PrintCerr(self,text):
        """To print to the original cerr"""
        self.stderr_.write(text+"\n")

def CheckIntegrity():
    #carefull, this class is used during the test.
    #do not use this class inside a CheckIntegrity
    from BasicTools.Helpers.Tests import TestTempDir
    fname = TestTempDir.GetTempPath() + "sink"
    with PrintBypass() as f:
        print("print Before ToSink")
        f.ToSink();
        f.Print("Print to the original cout")
        f.PrintCerr("Print to the original cerr")
        f.Restore()
        print("print after ToSink and befor ToFile")
        print("")

        f.ToDisk(filename=fname+"_cout.txt" )
        f.ToDisk(filename=fname+"_cout.txt",filenamecerr= fname+"_cerr.txt" )
        print("print inside to file cout")
        eprint("print inside to file cerr")

        f.Restore()
        print("print after ToFile and before custom class")

        print("")
        class mySink():
            def flush(self):pass
            def close(self): pass
            def write(self,data): pass

        f.ToRedirect(mySink(), mySink())
        eprint(' sent to cerr (mySink)')
        print('sent to cout (mySink)')

        f.Restore()
        print("print after mySink and before CopyOutputToDisk")
        print("")

        f.CopyOutputToDisk(filename=fname+"_cout_copy.txt" )
        f.CopyOutputToDisk(filename=fname+"_cout_copy.txt",filenamecerr= fname+"_cerr_copy.txt" )
        eprint(' sent to cerr (copy to disk)')
        print(' sent to cout (copy to disk) ')
        eprint(' sent to cerr (copy to disk) 2')
        print(' sent to cout (copy to disk) 2')
        f.Restore()



    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover