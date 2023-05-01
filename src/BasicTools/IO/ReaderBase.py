# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

"""Base reader object from which all the readers of BasicTools inherit
"""

import os
import sys
import struct
import locale

import numpy as np

from BasicTools.Helpers.BaseOutputObject import BaseOutputObject

class ReaderBase(BaseOutputObject):
    """ ReaderBase class"""
    def __init__(self,fileName = None)    :
        super(ReaderBase,self).__init__()
        self.fileName = None
        self.readFormat = 'r'
        self.nodalFields = {}
        self.elementsFields = {}
        self.binary = False
        self.string = None
        self.commentChar = None
        self.filePointer = None
        self.lineCounter = 0
        self.pipe = False
        self.encoding = locale.getpreferredencoding(False)
        self.canHandleTemporal = False

        self.output = None
        self.extraOutput = None

        self.SetFileName(fileName)

    def SetBinary(self, binary = True):
        """Sets the binary status of the file to read

        Parameters
        ----------
        binary : bool, optional
            if True, sets the file to read as binary, by default True
        """
        self.binary = binary
        if binary:
            if self.readFormat.find("b") >= 0:
                return
            self.readFormat += "b"
        else:
            if self.readFormat.find("b") >= 0:
                self.readFormat = self.readFormat.replace("b","")


    def StartReading(self):
        if not(self.fileName is None):
            if self.readFormat.find('b') > -1 :
                self.filePointer = open(self.fileName, self.readFormat)
                self.text_stream = self.filePointer
            else:
                #I have some problems reading with numpy fromfile if the file is
                #open with the codecs.open
                #import sys
                self.filePointer = open(self.fileName, self.readFormat,encoding=self.encoding)
                #import codecs
                #if sys.version_info[0] == 2:
                #    self.filePointer = codecs.open(self.fileName, self.readFormat, 'utf-8')
        elif not(self.string is None):
            if self.readFormat.find('b') > -1 :

                from io import BytesIO
                if type(self.string) is str:
                    self.filePointer = BytesIO(bytearray(self.string,self.encoding))
                else:
                    self.filePointer = BytesIO(self.string)
                self.text_stream = self.filePointer
            else:
                import io # Python3
                self.filePointer = io.StringIO(self.string)

                self.text_stream = self.filePointer
        elif self.pipe:
            import os
            r, w = os.pipe()
            if self.readFormat.find('b') > -1 :
                self.filePointer = sys.stdin.buffer
                #os.fdopen(r, self.readFormat)
                self.text_stream = self.filePointer
            else:
                #I have some problems reading with numpy fromfile if the file is
                #open with the codecs.open
                #import sys
                self.filePointer =  sys.stdin
                #os.fdopen(r, self.readFormat)
        else:
            raise ('Need a file or a string to read')

        self.lineCounter = 0

    def GetFilePointer(self):
        return self.filePointer

    def EndReading(self):
        self.filePointer.close()

    def SetFileName(self, fileName):
        """Sets the name of file to read

        Parameters
        ----------
        fileName : str
            file name to set
        """
        if  not(fileName is None) and len(fileName) >= 4 and fileName[0:4] == "PIPE" :
            self.SetReadFromPipe()
        else:
            self.fileName = fileName;
            if fileName is None :
                self.__path = None;
                self.string = None;
            else:
                self.filePath = os.path.abspath(os.path.dirname(fileName))+os.sep;
                self.string = None
                self.pipe = False

    def SetStringToRead(self,string):
        """Sets data to be read as a string instead of a file

        Parameters
        ----------
        string : str
            data to be read
        """
        self.string = string
        if string is not None:
            self.fileName = None
            self.pipe = False

    def SetReadFromPipe(self):
        self.SetFileName(None)
        self.SetStringToRead(None)
        self.pipe = True

    def PeekLine(self,withError=False):
        """Read a line without advancing"""

        pos = self.filePointer.tell()
        line = self.filePointer.readline()
        self.filePointer.seek(pos)
        return line

    def Peek(self,length=1):
        """Read a length number of chars without advancing the file """
        pos = self.filePointer.tell()
        data = self.filePointer.read(length) # Might try/except this line, and finally: f.seek(pos)
        self.filePointer.seek(pos)
        return data

    def ReadCleanLine(self,withError=False):
        while(True):
            string = self.filePointer.readline()

            self.lineCounter +=1
            #end of file
            if string == "" :
                if withError :
                    if self.fileName is None:
                        raise("Problem reading string : at line " +str(self.lineCounter))
                    else:
                        raise(Exception("Problem reading file :" +str(self.fileName) + " at line" +str(self.lineCounter) ))

                return None

            string = string.replace(u'\ufeff', '').strip(u' \r\n')

            #empty line
            if len(string) == 0 :
                continue
            if self.commentChar is None:
                break# pragma: no cover
            else :
                inbreak = False
                for i,j in zip(string,self.commentChar):
                    if i != j :
                        inbreak = True
                        break
                if inbreak:
                    break
                #if string[0] != self.commentChar:
                #    break
        return string

##binary interface
    def rawread(self,cpt,withError=False):

        res = self.filePointer.read(cpt)
        if withError and len(res) == 0:
            raise EOFError("Problem reading file :" +str(self.fileName) + " EOF")
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
            return np.frombuffer(data,dtype=datatype)

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
    obj.commentChar = "#"
    obj.SetBinary(False)
    try:
        obj.StartReading()
        raise # pragma: no cover
    except :
        pass
    testString = """0
1

2
3
#this is a comment
4"""
    obj.SetStringToRead(testString)

    def checkBaseReaderAscii(obj):
        obj.StartReading()
        print("file Pointer: ", str(obj.GetFilePointer() ) )
        if obj.PeekLine() != "0\n":
            raise
        if obj.Peek() != "0":
            raise

        for i in range(5):
            if i != int(obj.ReadCleanLine()):
                raise
        #before last
        if obj.ReadCleanLine() != None:
            raise

        error = False
        try:
            obj.ReadCleanLine(True) # this line must fail
            error = True
        except:
            pass

        if error:
            raise(Exception("Must fail ") )

        obj.EndReading()

    checkBaseReaderAscii(obj)
    from BasicTools.Helpers.Tests import WriteTempFile
    fn = WriteTempFile('ReaderBaseTestString',testString)
    obj.SetFileName(fn)
    checkBaseReaderAscii(obj)


    binarydata = np.array([0], dtype=np.int32).tobytes()
    binarydata += np.array([1], dtype=np.int64).tobytes()
    binarydata += np.array([2], dtype=np.float32).tobytes()
    binarydata += np.array([3], dtype=np.float64).tobytes()

    fn = WriteTempFile('ReaderBaseTestbinary',binarydata,mode='wb')
    obj.SetBinary(False)
    obj.SetBinary(True)
    obj.SetBinary(False)
    obj.SetBinary(False)
    obj.SetBinary(True)
    obj.SetBinary(True)
    obj.SetFileName(fn)
    print(obj.readFormat)
    def checkBaseReaderBinary(obj):
        obj.StartReading()
        if obj.readInt32() != 0: raise
        if obj.readInt64() != 1: raise
        if obj.readFloats32(1) != 2.: raise
        if obj.readFloats64(1) != 3.: raise
        obj.EndReading()
    checkBaseReaderBinary(obj)
    obj.SetStringToRead(binarydata)
    checkBaseReaderBinary(obj)

    return "ok"


if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
