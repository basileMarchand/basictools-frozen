# -*- coding: utf-8 -*-

import sys
import os
"""
Send all the print statements to a file or to the sink
"""

class PrintBypass():

    def __init__(self):
        self.stdout_ = sys.stdout #Keep track of the previous value.
        self.stderr_ = sys.stderr #Keep track of the previous value.
        self.bypassCout = False;
        self.bypassCerr = False;
        self.fileno = sys.stdout

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.Restore()

    def ToSink(self):
        sys.stdout = open(os.devnull,"w")
        sys.stderr = open(os.devnull,"w")
        self.bypassCout = True;
        self.bypassCerr = True;

    def ToDisk(self,filename, filenamecerr=None):
        sys.stdout = open(filename, 'w') # Something here that provides a write method.
        self.bypassCout = True;
        if filenamecerr is not  None:
            sys.stderr = sys.stdout
        else:
            sys.stderr = open(filenamecerr, 'w')
            self.bypassCerr = True;

    def ToRedirect(self,cout_obj,cerr_obj=None):
        # obj must implement the functions: close(), flush(), write(data)
        sys.stdout = cout_obj
        self.bypassCout = True;
        if cerr_obj is not None:
            sys.stderr = cerr_obj
            self.bypassCerr = True;

    def Restore(self):
        self.Print("Restore pipes")
        if self.bypassCout :
            sys.stdout.close()
        if self.bypassCerr :
            sys.stderr.close()
        sys.stdout = self.stdout_ # restore the previous stdout.
        sys.stderr = self.stderr_ # restore the previous stdout.

    def Print(self,text):
        """To print to the original cout"""
        self.stdout_.write(text+"\n")

    def PrintCerr(self,text):
        """To print to the original cerr"""
        self.stderr_.write(text+"\n")

def CheckIntegrity():
    #carefull, this class is used during the test.
    #do not use this class inside a CheckIntegrity
    PB = PrintBypass()
    return "ok"
