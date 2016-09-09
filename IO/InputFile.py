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
        self.filePointer = open(self.fileName,'r')

    def Close(self):
        self.filePointer.close()

    def GetCObjectName(self):
        pass

def CheckIntegrity():

    data = """ this is an input file """

    from OTTools.Helpers.Tests import TestTempDir
    tempdir = TestTempDir.GetTempPath()
    f =open(tempdir+"test_input_data","w")
    f.write(data)
    f.close()


    myObj = InputFile()
    myObj.SetFileName(tempdir+"test_input_data")
    myObj.Open()
    myObj.Close()
    myObj.GetCObjectName()

    return 'ok'

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
