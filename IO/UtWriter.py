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
        self.Nint             = data_integ[self.data_integ_names[0]].shape[0]

    def AttachSequence(self, cycle_number, sequence_number, increment, time):
        self.cycle_number = cycle_number
        self.sequence_number = sequence_number
        self.increment = increment
        self.time = time
        self.Ntime = len(time)

    def WriteMesh(self):
        if self.mesh==None:
          print("please attach a mesh to the UtWriter object to be able to write a mesh; script terminated")
          exit()
        else:
          OW = GW.GeofWriter()
          OW.Open(self.folder+self.name+".geof")
          OW.Write(self.mesh, useOriginalId=True)
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
        for i in range(len(self.increment)):
            __string += str(i+1)+" "+str(self.cycle_number[i])+" "+str(self.sequence_number[i])+" "+str(self.increment[i])+" "+str(self.time[i])+"\n"

        with open(self.folder+self.name+".ut", "w") as f:
          f.write(__string)
        f.close()


def CheckIntegrity():

    from BasicTools.Helpers.Tests import TestTempDir
    tempdir = TestTempDir.GetTempPath()
    import BasicTools.TestData as BasicToolsTestData

    cycle_number = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]
    sequence_number = [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
    increment = [0, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    time = [0., 0.2, 0.4, 0.6, 0.8, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.2, 2.4, 2.6, 2.8, 3., 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4., 4.2, 4.4, 4.6, 4.8, 5., 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.]

    import BasicTools.IO.UtReader as UR
    reader = UR.UtReader()
    reader.SetFileName(BasicToolsTestData.GetTestDataPath() + "UtExample/cube.ut")
    reader.ReadMetaData()

    reader.atIntegrationPoints = False
    Nnode = reader.Read(fieldname="U1", timeIndex=0).shape[0]
    reader.atIntegrationPoints = True
    Nint = reader.Read(fieldname="sig11", timeIndex=0).shape[0]
    Ntime = len(time)

    data_node_names = ['U1', 'U2', 'U3']
    data_integ_names = ['eto11', 'eto22', 'eto33', 'eto12', 'eto23', 'eto31', 'sig11', 'sig22', 'sig33', 'sig12', 'sig23', 'sig31']

    NnodeVar = len(data_node_names)
    NintVar = len(data_integ_names)

    import collections
    data_node = collections.OrderedDict()
    for dnn in data_node_names:
      data_node[dnn] = np.empty((Nnode,Ntime))
    data_ctnod = collections.OrderedDict()
    data_integ = collections.OrderedDict()
    for din in data_integ_names:
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
    UtW.AttachSequence(cycle_number, sequence_number, increment, time)
    UtW.Write(writeGeof=True)
    ##################################
    
    print(UtW)

    print("Temp directory =", tempdir)

    import filecmp
    print("node files equals  ?", filecmp.cmp(tempdir + "toto.node",  BasicToolsTestData.GetTestDataPath() + "UtExample/toto.node", shallow=False))
    print("ctnod files equals ?", filecmp.cmp(tempdir + "toto.ctnod", BasicToolsTestData.GetTestDataPath() + "UtExample/toto.ctnod", shallow=False))
    print("integ files equals ?", filecmp.cmp(tempdir + "toto.integ", BasicToolsTestData.GetTestDataPath() + "UtExample/toto.integ", shallow=False))
    print("ut files equals    ?", filecmp.cmp(tempdir + "toto.ut", BasicToolsTestData.GetTestDataPath() + "UtExample/toto.ut", shallow=False))

    return "ok"

if __name__ == '__main__':
    print((CheckIntegrity()))# pragma: no cover
