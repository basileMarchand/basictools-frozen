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
        self.name = name;

    def SetFolder(self,folder):
        self.folder = folder;

    def AttachMesh(self,mesh):
        self.mesh = mesh;

    def AttachData(self, data_node, data_ctnod, data_integ, data_node_names, data_integ_names):
        #data_node and data_integ are lists of numpy arrays, see example below
        self.data_node = data_node
        self.data_node_names = data_node_names
        self.data_ctnod = data_ctnod
        self.data_integ = data_integ
        self.data_integ_names = data_integ_names

    def AttachSequence(self, cycle_number, sequence_number, increment, time):
        self.cycle_number = cycle_number
        self.sequence_number = sequence_number
        self.increment = increment
        self.time = time

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

        try:
            self.data_node.astype(np.float32).byteswap().tofile(self.folder+self.name+".node")
            __string += "**node "
            for field in self.data_node_names:
              __string += field+" "
            __string += "\n"
        except AttributeError:
            print("non node data")

        try:
            self.data_ctnod.astype(np.float32).byteswap().tofile(self.folder+self.name+".ctnod")
        except AttributeError:
            print("non ctnod data")

        try:
            self.data_integ.astype(np.float32).byteswap().tofile(self.folder+self.name+".integ")
            __string += "**integ "
            for field in self.data_integ_names:
              __string += field+" "
            __string += "\n"
        except AttributeError:
            print("non integ data")

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

    data_node = np.empty(NnodeVar*Nnode*Ntime)
    data_ctnod = np.empty(NintVar*Nnode*Ntime)
    data_integ = np.empty(NintVar*Nint*Ntime)

    # CONVENTIONS POUR CTNOD ET NODE
    reader.atIntegrationPoints = False
    for i in range(Ntime):
      for k in range(NnodeVar):
        data_node[NnodeVar*Nnode*i+k*Nnode:NnodeVar*Nnode*i+(k+1)*Nnode] = reader.Read(fieldname=data_node_names[k], timeIndex=i)
      for k in range(NintVar):
        data_ctnod[NintVar*Nnode*i+k*Nnode:NintVar*Nnode*i+(k+1)*Nnode] = reader.Read(fieldname=data_integ_names[k], timeIndex=i)


    import BasicTools.IO.GeofReader as GR
    import BasicTools.FE.UnstructuredMesh as UM
    mymesh = GR.ReadGeof(fileName=BasicToolsTestData.GetTestDataPath() + "UtExample/cube.geof",out = UM.UnstructuredMesh())

    import BasicTools.FE.UnstructuredMeshTools as UnstructuredMeshTools
    from BasicTools.IO.GeofWriter import GeofName as GeofName
    from BasicTools.IO.GeofReader import nbIntegrationsPoints as nbIntegrationsPoints
    mymesh = UnstructuredMeshTools.ExtractElementByDimensionalityNoCopy(mymesh,3)
    numberElements = []
    nbPtIntPerElement = []
    for name,data in mymesh.elements.items():
      numberElements.append(data.GetNumberOfElements())
      nbPtIntPerElement.append(nbIntegrationsPoints[GeofName[name]])
    nbTypeEl = len(numberElements)



    # CONVENTIONS POUR .INTEG
    reader.atIntegrationPoints = True
    count0 = 0
    for i in range(Ntime):
      field = []
      for k in range(NintVar):
        field.append(reader.Read(fieldname=data_integ_names[k], timeIndex=i))
      for l in range(nbTypeEl):
        for m in range(numberElements[l]):
          for k in range(NintVar):
            data_integ[count0:count0+nbPtIntPerElement[l]] = field[k][nbPtIntPerElement[l]*m:nbPtIntPerElement[l]*m+nbPtIntPerElement[l]]
            count0 += nbPtIntPerElement[l]

    ##################################
    # EXEMPLE SYNTAXE DU WRITER  
    import BasicTools.IO.UtWriter as UW
    UtW = UW.UtWriter()
    UtW.SetName("toto")
    UtW.SetFolder(tempdir)
    UtW.AttachMesh(mymesh)
    UtW.AttachData(data_node, data_ctnod, data_integ, data_node_names, data_integ_names)
    UtW.AttachSequence(cycle_number, sequence_number, increment, time)
    UtW.Write(writeGeof=True)
    ##################################
    
    print("UtW =",UtW)

    print("tempdir =", tempdir)

    import filecmp
    print("node files equals  ?", filecmp.cmp(tempdir + "toto.node",  BasicToolsTestData.GetTestDataPath() + "UtExample/toto.node", shallow=False))
    print("ctnod files equals ?", filecmp.cmp(tempdir + "toto.ctnod", BasicToolsTestData.GetTestDataPath() + "UtExample/toto.ctnod", shallow=False))
    print("integ files equals ?", filecmp.cmp(tempdir + "toto.integ", BasicToolsTestData.GetTestDataPath() + "UtExample/toto.integ", shallow=False))
    print("ut files equals    ?", filecmp.cmp(tempdir + "toto.ut", BasicToolsTestData.GetTestDataPath() + "UtExample/toto.ut", shallow=False))

    return "ok"

if __name__ == '__main__':
    print((CheckIntegrity()))# pragma: no cover
