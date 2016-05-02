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
        self.bypass = False;
        self.fileno = sys.stdout
        
    def ToSink(self):
        sys.stdout = open(os.devnull,"w")
        sys.stderr = open(os.devnull,"w")
        self.bypass = True;

    def ToDisk(self,filename, filenamecerr=None):
        sys.stdout = open(filename, 'w') # Something here that provides a write method.    
        if filenamecerr is not  None:
            sys.stderr = sys.stdout
        else:
            sys.stderr = open(filenamecerr, 'w') # Something here that provides a write method.
        self.bypass = True;

    def Restore(self):
        if self.bypass :
            sys.stdout.close()
            sys.stderr.close()
        sys.stdout = self.stdout_ # restore the previous stdout.
        sys.stderr = self.stderr_ # restore the previous stdout.

    def Print(self,text):
        self.stdout_.write(text+"\n")

    def PrintCerr(self,text):
        self.stderr_.write(text+"\n")

def CheckIntegrity():
    #carefull, this class is used during the test.
    #do not use this class inside a CheckIntegrity
    PB = PrintBypass()
    return "ok"
