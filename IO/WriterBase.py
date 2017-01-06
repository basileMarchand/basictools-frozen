# -*- coding: utf-8 -*-
from OTTools.Helpers.BaseOutputObject import BaseOutputObject
from OTTools.Helpers.TextFormatHelper import TFormat as TFormat

class WriterBase(BaseOutputObject):
    def __init__(self, fileName = None):
        super(WriterBase,self).__init__()
        self.fileName = None;
        self.SetFileName(fileName)
        self._isOpen = False
        self._isBinary = False

    def SetBinary(self, val = True):
        if self._isOpen :
            print(TFormat.InRed("SetBinary before opening"))
            raise Exception
        self._isBinary = val

    def isBinary(self):
        return self._isBinary

    def isOpen(self):
        return self._isOpen


    def SetFileName(self,fileName):
        self.fileName = fileName;

    def Open(self, filename = None):
        if self._isOpen :
            print(TFormat.InRed("The file is already open !!!!!"))# pragma: no cover
            raise Exception


        if filename is not None:
            self.SetFileName(filename)

        ## we use unbuffered so we can repaire broken files easily
        try :
            if self._isBinary == True :
                mode = "wb"
            else:
                mode = "w"
            self.filePointer = open(self.fileName, mode,0)

        except:
            print(TFormat.InRed("Error File Not Open"))# pragma: no cover
            raise

        self._isOpen = True

    def Close(self):
        if self._isOpen:
            self.filePointer.close()
            self._isOpen = False
        else :
            self.PrintVerbose(TFormat.InRed("File Not Open"))

def CheckIntegrity():

    obj = WriterBase()
    return "ok"


if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
