# -*- coding: utf-8 -*-

from OTTools.Helpers.BaseOutputObject import BaseOutputObject

class InputFile(BaseOutputObject):
    def __init__(self):
        super(InputFile,self).__init__()
        self.fileName = None
        self.filePointer = None
        
    def SetFileName(self, fname):
        self.fileName = fname

    def Open(self):
        self.filePointer = open(self.filename,'r')
        
    def Close(self):
        self.filePointer.close()
        
    def GetCObjectName()
        
def CheckIntegrity():
    myObj = InputFile()
    
    
        
    return 'ok'
        
        
if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover 