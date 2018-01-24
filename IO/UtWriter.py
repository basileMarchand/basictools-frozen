# -*- coding: utf-8 -*-

#for python 2.6+ compatibility
from __future__ import print_function

import numpy as np
__author__ = "Fabien Casenave"

from BasicTools.IO.WriterBase import WriterBase as WriterBase
import BasicTools.IO.GeofWriter as GW


class UtWriter(WriterBase):
    "This class can generate a .ut, .goef, .ctnod, .node, .integ to vizualise finite element computational results"
    def __init__(self):
        super(UtWriter,self).__init__()
        self.canHandleTemporal = True
        self.canHandleAppend = True

    def __str__(self):
        res  = 'UtWriter : \n'
        res += '   Name : '+str(self.name)+'\n'
        return res

    def SetName(self,name):
        self.name = name

    def SetFolder(self,folder):
        self.folder = folder

    def AttachMesh(self,mesh):
        self.mesh = mesh

    def AttachData(self, data_node, data_ctnod = None, data_integ = None):
        self.data_node        = data_node
        self.data_ctnod       = data_ctnod
        self.data_integ       = data_integ
        self.data_node_names  = list(data_node.keys())
        self.data_integ_names = list(data_integ.keys())
        self.NnodeVar         = len(data_node)
        self.NintVar          = len(data_integ)
        self.Nnode            = data_node[self.data_node_names[0]].shape[0]
        try:
          self.Nint = data_integ[self.data_integ_names[0]].shape[0]
        except IndexError:
          self.Nint = None

    def AttachDataFromProblemData(self, problemData, tag):
        self.AttachData(problemData.solutions[tag].data_node, problemData.solutions[tag].data_ctnod, problemData.solutions[tag].data_integ)

    def AttachSequence(self, timeSequence):
        self.timeSequence = timeSequence
        self.Ntime = timeSequence.shape[0]

    def WriteMesh(self):
        if self.mesh==None:
          print("please attach a mesh to the UtWriter object to be able to write a mesh; script terminated")
          exit()
        else:
          OW = GW.GeofWriter()
          OW.Open(self.folder+self.name+".geof")
          OW.Write(self.mesh, useOriginalId=True, lowerDimElementsAsSets=True)
          OW.Close()

    def Write(self, writeGeof):

        __string = u"**meshfile "+self.name+".geof\n"

        if self.mesh==None:
          print("please attach a mesh to the UtWriter object to be able to write a mesh; script terminated")
          exit()

        if writeGeof==True:
          self.WriteMesh()

        if self.NnodeVar>0:
          data_node = np.empty(self.NnodeVar*self.Nnode*self.Ntime)
          for i in range(self.Ntime):
            for k in range(self.NnodeVar):
              data_node[self.NnodeVar*self.Nnode*i+k*self.Nnode:self.NnodeVar*self.Nnode*i+(k+1)*self.Nnode] = self.data_node[self.data_node_names[k]][:,i]
          data_node.astype(np.float32).byteswap().tofile(self.folder+self.name+".node")
          del data_node
          __string += "**node "
          for field in self.data_node_names:
            __string += field+" "
          __string += "\n"


        if self.NintVar>0:
          data_ctnod = np.empty(self.NintVar*self.Nnode*self.Ntime)
          for i in range(self.Ntime):
            for k in range(self.NintVar):
              data_ctnod[self.NintVar*self.Nnode*i+k*self.Nnode:self.NintVar*self.Nnode*i+(k+1)*self.Nnode] = self.data_ctnod[self.data_integ_names[k]][:,i]
          data_ctnod.astype(np.float32).byteswap().tofile(self.folder+self.name+".ctnod")
          del data_ctnod

          from BasicTools.IO.GeofWriter import GeofName as GeofName
          from BasicTools.IO.GeofReader import nbIntegrationsPoints as nbIntegrationsPoints
          import BasicTools.FE.UnstructuredMeshTools as UnstructuredMeshTools
          numberElements = []
          nbPtIntPerElement = []
          mesh3D = UnstructuredMeshTools.ExtractElementByDimensionalityNoCopy(self.mesh,3)
          for name,data in mesh3D.elements.items():
            numberElements.append(data.GetNumberOfElements())
            nbPtIntPerElement.append(nbIntegrationsPoints[GeofName[name]])
          nbTypeEl = len(numberElements)

          data_integ = np.empty(self.NintVar*self.Nint*self.Ntime)
          count0 = 0
          for i in range(self.Ntime):
            field = np.empty((self.NintVar,self.Nint))
            for k in range(self.NintVar):
              field[k,:] = self.data_integ[self.data_integ_names[k]][:,i]
            for l in range(nbTypeEl):
              for m in range(numberElements[l]):
                for k in range(self.NintVar):
                  data_integ[count0:count0+nbPtIntPerElement[l]] = field[k,nbPtIntPerElement[l]*m:nbPtIntPerElement[l]*m+nbPtIntPerElement[l]]
                  count0 += nbPtIntPerElement[l]
          data_integ.astype(np.float32).byteswap().tofile(self.folder+self.name+".integ")
          del field; del data_integ  
          __string += "**integ "
          for field in self.data_integ_names:
            __string += field+" "
          __string += "\n"

        __string += "**element\n"
        for i in range(self.Ntime):    
            __string += str(int(self.timeSequence[i,0]))+" "+str(int(self.timeSequence[i,1]))+" "+str(int(self.timeSequence[i,2]))+" "+str(int(self.timeSequence[i,3]))+" "+str(self.timeSequence[i,4])+"\n"

        with open(self.folder+self.name+".ut", "w") as f:
          f.write(__string)
        f.close()


def CheckIntegrity():

    from BasicTools.Helpers.Tests import TestTempDir
    tempdir = TestTempDir.GetTempPath()
    import BasicTools.TestData as BasicToolsTestData

    import BasicTools.IO.UtReader as UR
    reader = UR.UtReader()
    reader.SetFileName(BasicToolsTestData.GetTestDataPath() + "UtExample/cube.ut")
    reader.ReadMetaData()

    timeSequence = reader.time

    reader.atIntegrationPoints = False
    Nnode = reader.Read(fieldname="U1", timeIndex=0).shape[0]
    reader.atIntegrationPoints = True
    Nint = reader.Read(fieldname="sig11", timeIndex=0).shape[0]

    Ntime = reader.time.shape[0]
    NnodeVar = len(reader.node)
    NintVar = len(reader.integ)

    import collections
    data_node = collections.OrderedDict()
    for dnn in reader.node:
      data_node[dnn] = np.empty((Nnode,Ntime))
    data_ctnod = collections.OrderedDict()
    data_integ = collections.OrderedDict()
    for din in reader.integ:
      data_ctnod[din] = np.empty((Nnode,Ntime))
      data_integ[din] = np.empty((Nint,Ntime))

    reader.atIntegrationPoints = False
    for i in range(Ntime):
      for dn in data_node:
        data_node[dn][:,i] = reader.Read(fieldname=dn, timeIndex=i)
      for dc in data_ctnod:
        data_ctnod[dc][:,i] = reader.Read(fieldname=dc, timeIndex=i)
    reader.atIntegrationPoints = True
    for i in range(Ntime):
      for di in data_integ:
        data_integ[di][:,i] = reader.Read(fieldname=di, timeIndex=i)

    import BasicTools.IO.GeofReader as GR
    mymesh = GR.ReadGeof(fileName=BasicToolsTestData.GetTestDataPath() + "UtExample/cube.geof")
  
    ##################################
    # EXEMPLE SYNTAXE DU WRITER  
    import BasicTools.IO.UtWriter as UW
    UtW = UW.UtWriter()
    UtW.SetName("toto")
    UtW.SetFolder(tempdir)
    UtW.AttachMesh(mymesh)
    UtW.AttachData(data_node, data_ctnod, data_integ)
    UtW.AttachSequence(timeSequence)
    UtW.Write(writeGeof=True)
    ##################################
    
    print(UtW)

    print("Temp directory =", tempdir)

    import filecmp
    print("node files equals  ?", filecmp.cmp(tempdir + "toto.node",  BasicToolsTestData.GetTestDataPath() + "UtExample/cube.node", shallow=False))
    print("ctnod files equals ?", filecmp.cmp(tempdir + "toto.ctnod", BasicToolsTestData.GetTestDataPath() + "UtExample/cube.ctnod", shallow=False))
    print("integ files equals ?", filecmp.cmp(tempdir + "toto.integ", BasicToolsTestData.GetTestDataPath() + "UtExample/cube.integ", shallow=False))

    return "ok"

if __name__ == '__main__':
    print((CheckIntegrity()))# pragma: no cover
