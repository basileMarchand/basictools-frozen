# -*- coding: utf-8 -*-
import os

from OTTools.Helpers.BaseOutputObject import BaseOutputObject
#from OTTools.Helpers.TextFormatHelper import TFormat as TFormat


class ReaderBase(BaseOutputObject):
    def __init__(self,fileName = None)    :
        super(ReaderBase,self).__init__()
        self.fileName = None;
        self.SetFileName(fileName)
        self.nodalFields = {}
        self.elementsFields = {}
        self.readFormat = 'r'
        self.output = None
        self.string = None

    def StartReading(self):

        if self.fileName is None:
            if self.string is None:
                raise ('Need a file or a string to read')
            else:
                from cStringIO import StringIO
                self.filePointer =  StringIO(self.string)
        else:
            self.filePointer =  open(self.fileName, self.readFormat)



    def EndReading(self):
        self.filePointer.close()

    def SetFileName(self,fileName):


        self.fileName = fileName;
        if fileName is None :
            self.__path = None;
            self.string = None;
        else:
            self.filePath = os.path.abspath(os.path.dirname(fileName));
            self.string = None


    def SetStringToRead(self,string):
        self.string = string

        if string is not None:
            self.fileName = None

def CheckIntegrity():

    obj = ReaderBase()
    try:
        obj.StartReading()
        raise # pragma: no cover
    except :
        pass
    return "ok"


if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
