# -*- coding: utf-8 -*-
import numpy as np

from BasicTools.IO.ReaderBase import ReaderBase


def ReadUt(fileName=None,fieldname=None,time=None,string=None,atIntegrationPoints=False):
    reader = UtReader()
    reader.SetFileName(fileName)
    reader.SetStringToRead(string)
    reader.atIntegrationPoints = atIntegrationPoints
    reader.ReadMetaData()
    return reader.Read(fieldname=fieldname, time=time )

def ReadMeshAndUt(fileName):
        reader = UtReader()
        reader.SetFileName(fileName=fileName)
        reader.ReadMetaData()
        import BasicTools.IO.GeofReader as GeofReader
        mesh = GeofReader.ReadGeof(reader.meshfile)

        #mesh.nodeFields = {}
        nodesfields = reader.node
        nodesfields.extend(reader.integ)
        for nf in nodesfields:
            mesh.nodeFields[nf] = reader.Read(fieldname=nf,time=-1)
        return mesh

class UtReader(ReaderBase):
    def __init__(self):
        super(UtReader,self).__init__()
        self.commentChar ="#"
        self.meshfile = None
        self.node = None
        self.integ = None
        self.element =None
        self.time = None

        self.fieldNameToRead = None
        self.timeToRead = -1
        self.atIntegrationPoints = False
        self.meshMetadata = None
        self.cache = None
        self.oldtimeindex = None

        self.canHandleTemporal = True

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

    def GetAvilableTimes(self):
           return self.time[:,4]

    def ReadMetaData(self):
        if self.meshMetadata is not None : return self.meshMetadata
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
        self.time = np.array(self.time)
        self.EndReading()

        from BasicTools.IO.GeofReader import GeofReader
        if self.meshfile[-5:] == ".geof":
            GR = GeofReader()
            GR.SetFileName(self.filePath +self.meshfile )
            self.meshMetadata = GR.ReadMetaData()
        else:
            self.Print("Unable to obtain metadata from meshfile, please set metadata manually")


    def ReadBinaryFile(self,fileName):
        return np.fromfile(fileName, dtype=np.float32).byteswap()


    def ReadField(self,fieldname=None,time=None,timeIndex=None):
        self.ReadMetaData()
        postfix = ""
        if self.fileName[-1] == "p":
            postfix = "p"

        self.SetFieldNameToRead(fieldname)

        if time is None:
          if timeIndex is None:
            timeIndex = len(self.time)-1
          else:
            timeIndex = timeIndex
        else:
          self.SetTimeToRead(time)

          # find the time
          if self.timeToRead == -1 :
            timeIndex = len(self.time)-1
          else:
            timeIndex = [data[4]for data in self.time].index(self.timeToRead)

        if timeIndex != self.oldtimeindex:
            self.cache = None
            self.oldtimeindex = timeIndex


        self.PrintVerbose("Reading timeIndex : " + str(timeIndex) )

        basename = ".".join(self.fileName.split(".")[0:-1])
        #find the field to read
        idx = None
        res = None

        nbUsednodes = self.meshMetadata['nbNodes']
        nbNodes = self.meshMetadata['nbNodes']
        nbIntegrationPoints = self.meshMetadata['nbIntegrationPoints']
        nbElements = self.meshMetadata['nbElements']
        IPPerElement = self.meshMetadata['IPPerElement']

        if self.atIntegrationPoints :
            try:
                idx = self.integ.index(self.fieldNameToRead)
                offset =  nbIntegrationPoints * len(self.integ)*timeIndex
                count = nbIntegrationPoints
                ffn = basename + ".integ"
            except:
                raise(Exception("unable to find field " +str(self.fieldNameToRead) ))

            self.PrintVerbose("Opening file : " + str(ffn) )
            res = np.empty(count,dtype=np.float)
            try:
              if len(self.integ)==1 :
                with open(ffn+postfix,"rb") as datafile:
                   datafile.seek(offset*4)
                   res = np.fromfile(datafile ,count=count, dtype=np.float32).byteswap()
              elif np.min(IPPerElement) == np.max(IPPerElement)  :
                  # the .intef file is homogenius
                with open(ffn+postfix,"rb") as datafile:
                    datafile.seek(offset*4)
                    nip = IPPerElement[0]

                    if self.cache is None:
                        self.cache = np.fromfile(datafile ,count=nip*nbElements* len(self.integ), dtype=np.float32).byteswap()
                        self.cache.shape = (nbElements,nip* len(self.integ))

                    res = self.cache[:,idx*nip:(idx+1)*nip].flatten()
              else:
                with open(ffn+postfix,"rb") as datafile:
                    self.PrintDebug("Offset : " + str(offset*4))
                    self.PrintDebug("count : " + str(count))
                    datafile.seek(offset*4)
                    cpt =0
                    oldStep = 0
                    for el in range(nbElements):
                        #for sv in range(len(self.integ)):
                            #for ip in range(IPPerElement[el]):
                                #res[cpt] = np.fromfile(datafile ,count=1, dtype=np.float32).byteswap()

                        nip = IPPerElement[el]
                        oldStep += nip*idx
                        datafile.seek(4*(oldStep),1)
                        res[cpt:cpt +nip] = np.fromfile(datafile ,count=nip, dtype=np.float32).byteswap()
                        cpt += nip

                        oldStep = (len(self.integ)-idx-1)*nip

            except:
                print("Error Reading field : " + str(self.fieldNameToRead) + " (not read)")

        else:
            try:
                idx = self.node.index(self.fieldNameToRead)
                offset = nbNodes * len(self.node)*timeIndex + idx*nbNodes
                count = nbNodes
                ffn = basename + ".node"
            except:
                try:
                    idx = self.integ.index(self.fieldNameToRead)
                    offset = nbUsednodes * len(self.integ)*timeIndex + idx*nbUsednodes
                    count = nbUsednodes
                    ffn = basename + ".ctnod"
                except:
                    raise(Exception("unable to find field " +str(self.fieldNameToRead) + "in file " + self.fileName))


            self.PrintVerbose("Opening file : " + str(ffn) )
            res = None
            try:
                with open(ffn+postfix,"rb") as datafile:
                    self.PrintDebug("Offset : " + str(offset*4))
                    self.PrintDebug("count : " + str(count))
                    datafile.seek(offset*4)
                    res = np.fromfile(datafile ,count=count, dtype=np.float32).byteswap()
            except:
                print("Error Reading field : " + str(self.fieldNameToRead) + " (not read)")

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
    nbElements = reader.meshMetadata['nbElements']
    IPPerElement = reader.meshMetadata['IPPerElement']

    np.arange(nbNodes*6*2, dtype=np.float32).byteswap().tofile(TestTempDir.GetTempPath() + "UtReaderTest.node")
    off1 = nbNodes*6*2
    (np.arange(nbNodes*12*2, dtype=np.float32)+off1).byteswap().tofile(TestTempDir.GetTempPath() + "UtReaderTest.ctnod")
    off2 = (nbNodes*6*2+nbNodes*12*2)

    ipdata = np.empty((nbElements,12,27)  ,dtype=int)
    for el in range(nbElements):
         nip = IPPerElement[el]
         for sv in range(12):
            for ip in range(nip):
                     ipdata[el,sv,ip] = el*10000+ sv*100+ip
#                    res[cpt] = np.fromfile(datafile ,count=1, dtype=np.float32).byteswap()

    #print(ipdata[:,0,:] )
    #ipdata.astype(np.float32).byteswap().tofile(TestTempDir.GetTempPath() + "UtReaderTest.integ")

    #(np.arange(216*12*27*2, dtype=np.float32)+off2).byteswap().tofile(TestTempDir.GetTempPath() + "UtReaderTest.integ")


    with open(TestTempDir.GetTempPath() + "UtReaderTest.integ","wb") as datafile:
      for i in range(2):
        for el in range(nbElements):
            nip = IPPerElement[el]
            for sv in range(12):
                for ip in range(nip):
                    np.array(ipdata[el,sv,ip] + i *1000000).astype(np.float32).byteswap().tofile(datafile)
#                    res[cpt] = np.fromfile(datafile ,count=1, dtype=np.float32).byteswap()
#
#

    offset = 0
    for t in [0., 1.]:
        for f in ["U1","U2","U3","RU1","RU2","RU3"]:
            if np.any (reader.Read(fieldname=f,time=t) != np.arange(offset, offset+nbNodes, dtype=np.float32) ):
                raise(Exception("Error Reading field " + f ))
            offset += nbNodes

    offset = 0
    for t in [0., 1.]:
        for f in ["sig11","sig22","sig33","sig12","sig23","sig31","eto11","eto22","eto33","eto12","eto23","eto31"]:
            data = reader.Read(fieldname=f,time=t)
            if np.any (data !=  np.arange(offset, offset+nbNodes, dtype=np.float32)+off1 ):
                raise(Exception("Error Reading field " + f ))
            offset += nbNodes

    reader.atIntegrationPoints = True

    offset = 0
    for t in [0., 1.]:
        cpt =0
        for f in ["sig11","sig22","sig33","sig12","sig23","sig31","eto11","eto22","eto33","eto12","eto23","eto31"]:

            data = reader.Read(fieldname=f,time=t)
            if np.any (data != ipdata[:,cpt,:].ravel()+ t *1000000 ):
                print(data)
                print(' '.join( [ str(int(x)) for x in data[0:30] ] ))
                print(' '.join( [ str(x) for x in ipdata[:,cpt,:].ravel()[0:30] ] ))
                print(len(data)),
                print(ipdata[:,cpt,:].size),
                print(ipdata.shape)

                raise(Exception("Error Reading field " + f ))
            offset += nbIP
            cpt +=1


    #res = ReadUt(string=__teststring,fieldname="U1")

    ReadUt(tempfileName,fieldname="U1",time=0)

    return "OK"

if __name__ == '__main__':# pragma: no cover
    #from BasicTools.Helpers.BaseOutputObject import BaseOutputObject
    #BaseOutputObject.SetGlobalDebugMode(True)
    import time
    a =  time.time()
    print(CheckIntegrity())# pragma: no cover
    print(time.time() -a)
