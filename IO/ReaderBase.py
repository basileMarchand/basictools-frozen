# -*- coding: utf-8 -*-
import os
import struct

import numpy as np

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
        self.lineCounter = 0

    def StartReading(self):

        import sys

        if self.fileName is None:
            if self.string is None:
                raise ('Need a file or a string to read')
            else:
                if self.readFormat.find('b') > -1 :

                    from io import BytesIO
                    if sys.version_info >= (3,0):
                        self.filePointer =  BytesIO(bytearray(self.string,"ascii"))
                    else:
                        self.filePointer =  BytesIO(str(self.string))

                    self.text_stream = self.filePointer
                else:
                    #python2 test
                    if sys.version_info >= (3,0):
                        import io
                        self.filePointer =  io.StringIO(self.string)

                    else:
                        import cStringIO
                        self.filePointer =  cStringIO.StringIO(self.string)

                    self.text_stream = self.filePointer
        else:
            if self.readFormat.find('b') > -1 :
                self.filePointer =  open(self.fileName, self.readFormat)
                self.text_stream = self.filePointer
            else:
                #I have some problems reading with numpy fromfile if the file is
                #open with the codecs.open
                #import sys
                self.filePointer =  open(self.fileName, self.readFormat)
                #import codecs
                #if sys.version_info[0] == 2:
                #    self.filePointer = codecs.open(self.fileName, self.readFormat, 'utf-8')




        self.lineCounter = 0

    def GetFilePointer(self):
         return self.filePointer

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

    def ReadCleanLine(self,withError=False):
        while(True):
            import sys
            if sys.version_info[0] == 2:
                string = self.filePointer.readline().decode("utf-8", "replace")
            else:
                string = self.filePointer.readline()
            #print(string)
            #import codecs
            #string = codecs.getreader("utf-8")(self.filePointer).readline()
            #
            # old code working only in python 2

            self.lineCounter +=1
            #end of file
            if string == "" :
                if withError :
                    if self.fileName is None:
                        raise("Problem reading string : at line " +str(self.lineCounter))
                    else:
                        raise(Exception("Problem reading file :" +str(self.fileName) + " at line" +str(self.lineCounter) ))

                return None

            string = string.replace(u'\ufeff', '')
            string = string.lstrip().rstrip().rstrip(u'\r\n')
            #empty line
            if len(string) == 0 :
                continue
            if self.commentChar is None:
                break# pragma: no cover
            else :
                if string[0] != self.commentChar:
                    break
        return string

##binary interface
    def rawread(self,cpt,withError=False):

        res = self.filePointer.read(cpt)
        if withError and len(res) == 0:
           raise(Exception("Problem reading file :" +str(self.fileName) + " EOF"))
        else:
           return res

    def readInt32(self):
       rawdata = self.rawread(4,withError=True)
       data = struct.unpack("i", rawdata)[0]
       return data

    def readInt64(self):
       rawdata = self.rawread(8,withError=True)
       data = struct.unpack("q", rawdata)[0]
       return data

    def readData(self,cpt,datatype):
        try:
            return  np.fromfile(self.filePointer,dtype=datatype,count=cpt,sep="")
        except:
            s = np.dtype(datatype).itemsize*cpt
            data = self.filePointer.read(s)
            return np.fromstring(data,dtype=datatype)

    def reshapeData(self,data,finalShape=None):
        if finalShape is None:
            return data
        else:
            data.shape = finalShape
            return data

    def readFloats32(self,cpt,finalShape=None):
        return self.reshapeData(self.readData(cpt,np.float32), finalShape)

    def readFloats64(self,cpt,finalShape=None):
        return self.reshapeData(self.readData(cpt,np.float64), finalShape)

    def seek(self,cpt):
        self.filePointer.seek(cpt)

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
