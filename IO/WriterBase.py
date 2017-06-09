# -*- coding: utf-8 -*-
from BasicTools.Helpers.BaseOutputObject import BaseOutputObject
from BasicTools.Helpers.TextFormatHelper import TFormat as TFormat

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
        if self._isOpen :# pragma: no cover
            print(TFormat.InRed("The file is already open !!!!!"))
            raise Exception


        if filename is not None:
            self.SetFileName(filename)

        ## we use unbuffered so we can repaire broken files easily
        try :
            if self._isBinary  :
                mode = "wb"
            else:
                # in python 3 the binary mode must be used to use the numpy.savetxt
                mode = "w"

            # unbuffered text I/O are not allowed in python 3
            # bug http://bugs.python.org/issue17404
            #import io
            #import sys
            #binstdout = io.open(self.fileName, 'wb', 0)
            #self.filePointer = io.TextIOWrapper(binstdout, encoding=sys.stdout.encoding)
            #
            self.filePointer = open(self.fileName, mode)

        except:# pragma: no cover
            print(TFormat.InRed("Error File Not Open"))# pragma: no cover
            raise

        self._isOpen = True

    def writeText(self,text):
        if self._isBinary:
            self.filePointer.write(text.encode('utf8'))
        else:
            self.filePointer.write(text)

    def Close(self):
        if self._isOpen:
            self.filePointer.close()
            self._isOpen = False
        else :
            self.PrintVerbose(TFormat.InRed("File Not Open"))
            raise

def CheckIntegrity():

    obj = WriterBase()
    return "ok"


if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
