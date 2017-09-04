# -*- coding: utf-8 -*-
import os

__author__ = "Felipe Bordeu"
from BasicTools.Helpers.BaseOutputObject import BaseOutputObject

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
        self.commentChar = None
        self.filePointer = None

    def StartReading(self):

        if self.fileName is None:
            if self.string is None:
                raise ('Need a file or a string to read')
            else:
                from io import StringIO
                self.filePointer =  StringIO(self.string)
        else:
            if self.readFormat.find('b') > -1 :
                self.filePointer =  open(self.fileName, self.readFormat)
            else:
                import codecs
                self.filePointer = codecs.open(self.fileName, self.readFormat, 'utf-8')

    def EndReading(self):
        self.filePointer.close()

    def SetFileName(self,fileName):

        self.fileName = fileName;
        if fileName is None :
            self.__path = None;
            self.string = None;
        else:
            self.filePath = os.path.abspath(os.path.dirname(fileName))+os.sep;
            self.string = None


    def SetStringToRead(self,string):
        self.string = string

        if string is not None:
            self.fileName = None

    def ReadCleanLine(self):
        while(True):
            string = self.filePointer.readline()
            #end of file
            if string == "" :
                return None

            string = string.rstrip(u'\r\n')
            string = string.replace(u'\ufeff', '').lstrip().rstrip()
            #empty line
            if len(string) == 0 :
                continue
            if self.commentChar is None:
                break# pragma: no cover
            else :
                if string[0] != self.commentChar:
                    break
        return string

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
