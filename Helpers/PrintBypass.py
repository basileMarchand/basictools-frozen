# -*- coding: utf-8 -*-

import sys

""" 
Send all the print statements to a file or to the sink
"""

class PrintBypass():
    class Sink():
        def write(self,tosink):
            pass
        def close(self):
            pass
        
    def __init__(self):
        self.stdout_ = sys.stdout #Keep track of the previous value.
        self.bypass = False;
        
    def ToSink(self):
        sys.stdout = PrintBypass().Sink()
        self.bypass = True;
        
    def ToDisk(self,filename):
        sys.stdout = open(filename, 'w') # Something here that provides a write method.    
        self.bypass = True;
        
    def Restore(self):
        if self.bypass :
            sys.stdout.close()
        sys.stdout = self.stdout_ # restore the previous stdout.      
        
    def Print(self,text):
        self.stdout_.write(text+"\n")
        
        