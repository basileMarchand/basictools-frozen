# -*- coding: utf-8 -*-
import numpy as np

from BasicTools.IO.ReaderBase import ReaderBase


def ReadUt(fileName=None,fieldname=None,time=None,string=None):
    reader = UtReader()
    reader.SetFileName(fileName)
    reader.SetStringToRead(string)
    return reader.Read(fieldname=fieldname, time=time )

class UtReader(ReaderBase):
    def __init__(self):
        super(UtReader,self).__init__()
        self.commentChar ="#"
        self.meshfile = None
        self.node = None
        self.integ = None
        self.element =None
        self.time = None

        self.fieldNameToRead =None
        self.timeToRead = -1
        self.atIntegrationPoints = False

    def Reset(self):
        self.meshfile = None
        self.node = None
        self.integ = None
        self.element =None
        self.time = []

    def SetFieldNameToRead(self,fieldname):
        if fieldname is not None:
            self.fieldNameToRead = fieldname

    def SetTimeToRead(self,time):
        if time is not None:
            self.timeToRead = time

    def ReadMetaData(self):
        self.StartReading()

        self.Reset()


        while(True):
            line = self.ReadCleanLine()
            self.PrintVerbose(line)
            if line == None :
                break
            if line.find("**meshfile")>-1 :
                s = line.split()
                if len(s) == 1 :
                    line = self.ReadCleanLine()
                    self.meshfile = line
                else:
                    self.meshfile = s[-1]
                continue

            if line.find("**node")>-1 :
                s = line.split()
                self.node = s[1:]
                continue

            if line.find("**integ")>-1 :
                s = line.split()
                self.integ = s[1:]
                continue

            if line.find("**element")>-1 :
                s = line.split()
                self.element = s[1:]
                continue

            #all the rest are the times
            s = line.split()

            self.time.append( [b(a) for a,b in zip(s,[int, int, int, int, float] )])
        self.EndReading()

        from BasicTools.IO.GeofReader import GeofReader
        GR = GeofReader()
        GR.SetFileName(self.filePath +self.meshfile )
        self.meshMetadata = GR.ReadMetaData()

    def Read(self,fieldname=None,time=None):
        self.ReadMetaData()
        postfix = ""
        if self.fileName[-1] == "p":
            postfix = "p"

        self.SetFieldNameToRead(fieldname)
        self.SetTimeToRead(time)

        # find the time
        if self.timeToRead is -1 or self.timeToRead is -1.:
            timeIndex = len(self.time)-1
        else:
            timeIndex = [data[4]for data in self.time].index(self.timeToRead)
        self.PrintVerbose("Reading timeIndex : " + str(timeIndex) )

        basename = ".".join(self.fileName.split(".")[0:-1])
        #find the field to read
        ok = False
        idx = None
        res = None


        nbUsednodes = self.meshMetadata['nbNodes']
        nbNodes = self.meshMetadata['nbNodes']
        nbIntegrationPoints = self.meshMetadata['nbIntegrationPoints']

        if self.atIntegrationPoints :
            try:
                idx = self.integ.index(fieldname)
                offset =  nbIntegrationPoints * len(self.integ)*timeIndex + idx*nbIntegrationPoints
                count = nbIntegrationPoints
                ffn = basename + ".integ"
            except:
                raise(Exception("unable to find field " +str(fieldname) ))
        else:
            try:
                idx = self.node.index(fieldname)
                offset = nbNodes * len(self.node)*timeIndex + idx*nbNodes
                count = nbNodes
                ffn = basename + ".node"
            except:
                try:
                    idx = self.integ.index(fieldname)
                    offset = nbUsednodes * len(self.integ)*timeIndex + idx*nbUsednodes
                    count = nbUsednodes
                    ffn = basename + ".ctnod"
                except:
                    raise(Exception("unable to find field " +str(fieldname) ))


        self.PrintVerbose("Opening file : " + str(ffn) )
        with open(ffn+postfix,"rb") as datafile:
            self.PrintDebug("Offset : " + str(offset*4))
            self.PrintDebug("count : " + str(count))
            datafile.seek(offset*4)
            res = np.fromfile(datafile ,count=count, dtype=np.float32).byteswap()
        return res

    def __str__(self):
        res = ""
        res +=  "class UtReader ("+str(id(self)) + ")\n"
        res +=  "  meshfile : "+str(self.meshfile)+ "\n"
        res +=  "  node : "+str(self.node)+ "\n"
        res +=  "  integ : "+str(self.integ)+ "\n"
        res +=  "  element : "+str(self.element)+ "\n"
        res +=  "  times : \n"
        for i in self.time:
            res += "     " +str(i)+ "\n"
        return res

def CheckIntegrity():

    __teststring = u"""
**meshfile cube.geof
**node U1 U2 U3 RU1 RU2 RU3
**integ sig11 sig22 sig33 sig12 sig23 sig31 eto11 eto22 eto33 eto12 eto23 eto31
**element
1 1 1 1 0.000000000000000e+00
1 1 1 1 1.000000000000000e+00
"""

    from BasicTools.Helpers.Tests import WriteTempFile
    tempfileName = WriteTempFile("UtReaderTest.ut",content=__teststring)


    import BasicTools.TestData as BasicToolsTestData
    from BasicTools.Helpers.Tests import TestTempDir

    import shutil

    shutil.copy(BasicToolsTestData.GetTestDataPath() + "cube.geof",
                TestTempDir.GetTempPath() + "cube.geof")
    reader = UtReader()
    reader.SetFileName(tempfileName)
    reader.ReadMetaData()
    nbNodes = reader.meshMetadata['nbNodes']
    nbIP = reader.meshMetadata['nbIntegrationPoints']

    np.arange(nbNodes*6*2, dtype=np.float32).byteswap().tofile(TestTempDir.GetTempPath() + "UtReaderTest.node")
    off1 = nbNodes*6*2
    (np.arange(nbNodes*12*2, dtype=np.float32)+off1).byteswap().tofile(TestTempDir.GetTempPath() + "UtReaderTest.ctnod")
    off2 = (nbNodes*6*2+nbNodes*12*2)
    (np.arange(216*12*27*2, dtype=np.float32)+off2).byteswap().tofile(TestTempDir.GetTempPath() + "UtReaderTest.integ")

    offset = 0
    for t in [0., 1.]:
        for f in ["U1","U2","U3","RU1","RU2","RU3"]:
            if np.any (reader.Read(fieldname=f,time=t) != np.arange(offset, offset+nbNodes, dtype=np.float32) ):
                raise("Error Reading field" + f )
            offset += nbNodes

    offset = 0
    for t in [0., 1.]:
        for f in ["sig11","sig22","sig33","sig12","sig23","sig31","eto11","eto22","eto33","eto12","eto23","eto31"]:
            data = reader.Read(fieldname=f,time=t)
            if np.any (data != np.arange(offset, offset+nbNodes, dtype=np.float32)+off1 ):
                raise("Error Reading field" + f )
            offset += nbNodes


    reader.atIntegrationPoints = True

    offset = 0
    for t in [0., 1.]:
        for f in ["sig11","sig22","sig33","sig12","sig23","sig31","eto11","eto22","eto33","eto12","eto23","eto31"]:

            data = reader.Read(fieldname=f,time=t)
            if np.any (data != np.arange(offset, offset+nbIP, dtype=np.float32)+off2 ):
                raise("Error Reading field" + f )
            offset += nbIP


    #res = ReadUt(string=__teststring,fieldname="U1")

    ReadUt(tempfileName,fieldname="U1",time=0)

    return "OK"

if __name__ == '__main__':# pragma: no cover
    #from BasicTools.Helpers.BaseOutputObject import BaseOutputObject
    #BaseOutputObject.SetGlobalDebugMode(True)
    print(CheckIntegrity())# pragma: no cover