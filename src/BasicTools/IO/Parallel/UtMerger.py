# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#


import numpy as np
import os

from BasicTools.IO.WriterBase import WriterBase as WriterBase
import BasicTools.IO.GeofWriter as GW
import BasicTools.IO.GeofReader as GR
import BasicTools.IO.UtReader as UR
from BasicTools.IO.GeofWriter import GeofName as GeofName
from BasicTools.IO.GeofReader import nbIntegrationsPoints as nbIntegrationsPoints
import BasicTools.Containers.ElementNames as EN
from BasicTools.FE.IntegrationsRules import Lagrange as Lagrange
from BasicTools.Helpers.TextFormatHelper import TFormat

class UtMerger(WriterBase):
    "This class can generate monolithic .ut, .goef, .ctnod, .node, .integ from files obtained from a parallel computation"
    def __init__(self):
        super(UtMerger,self).__init__()
        self.timeSteps = "all"

    def __str__(self):
        res  = 'UtMerge : \n'
        res += '   Name : '+str(self.name)+'\n'
        return res

    def SetName(self,name):
        self.name = name

    def SetdataFolder(self,folder):
        self.dataFolder = folder
        from os import listdir
        from os.path import isfile, join
        temp = [f.split(".") for f in listdir(self.dataFolder) if isfile(join(self.dataFolder, f))]
        temp2 = []
        count = 0
        for f in temp:
          temp2.append([])
          for fi in f:
            splitted = fi.split("-")
            for fil in splitted:
              try:
                temp2[count].append(int(fil))
              except ValueError:
                temp2[count].append(fil)
          count += 1
        subdomains = []
        for f in temp2:
          if 'ut' in f and self.name in f:
            for fi in f:
              if type(fi) == int:
                subdomains.append(fi)
        self.nbsd = np.max(subdomains)
        print("Number of found subdomains:", self.nbsd)

    def SetOutputFolder(self,outputFolder):
        self.outputFolder = outputFolder

    def SetTimeSteps(self,iterator):
        self.timeSteps = iterator

    def Merge(self):
      localDataInteg = []
      localDataNode = []
      localIdstotreat = []
      localoriginalIDNodes = []

      #Read each subdomain computation
      for sd in range(1,self.nbsd+1):
        if sd<10:
          sdString = "-00"+str(sd)
        elif sd<100:
          sdString = "-0"+str(sd)
        else:
          sdString = "-"+str(sd)

        reader = UR.UtReader()
        reader.SetFileName(self.dataFolder + self.name + sdString + ".ut")
        reader.ReadMetaData()

        if self.timeSteps != "all":
          reader.time = reader.time[self.timeSteps,:]
          if len(reader.time.shape) == 1:
            reader.time.shape = (1,reader.time.shape[0])

        localMesh = GR.ReadGeof(fileName = self.dataFolder+reader.meshfile,readElset=False,readFaset=False,printNotRead=False)
        Tag3D(localMesh)
        idstotreat = Return3DElements(localMesh)

        originalIDNodes = np.array(localMesh.originalIDNodes-1, dtype=int)

        reader.atIntegrationPoints = True
        dataInteg = {}
        for din in reader.integ:
          dataInteg[din] = np.empty((reader.meshMetadata['nbIntegrationPoints'],reader.time.shape[0]))
          for timeStep in range(reader.time.shape[0]):
            dataInteg[din][:,timeStep] = reader.ReadField(fieldname=din, timeIndex=int(reader.time[timeStep,0])-1)

        reader.atIntegrationPoints = False
        dataNode = {}
        for din in reader.node:
          dataNode[din] = np.empty((reader.meshMetadata['nbNodes'],reader.time.shape[0]))
          for timeStep in range(reader.time.shape[0]):
            dataNode[din][:,timeStep] = reader.ReadField(fieldname=din, timeIndex=int(reader.time[timeStep,0])-1)

        localDataInteg.append(dataInteg)
        localDataNode.append(dataNode)
        localIdstotreat.append(idstotreat)
        localoriginalIDNodes.append(originalIDNodes)

      cutGeof = GeofFromCut(self.dataFolder, self.name)
      globalMesh = GR.ReadGeof(fileName = self.dataFolder + cutGeof,readElset=False,readFaset=False,printNotRead=False)

      Tag3D(globalMesh)
      globalIdstotreat = Return3DElements(globalMesh)
      #globalMesh.originalIDNodes = np.array(globalMesh.originalIDNodes-1, dtype=int)

      globaldataInteg = {}
      for din in reader.integ:
        globaldataInteg[din] = np.empty((globalMesh.NGauss,reader.time.shape[0]))
      globaldataNode = {}
      for din in reader.node:
        globaldataNode[din] = np.empty((globalMesh.Nodes,reader.time.shape[0]))

      #write .integ
      data_integ = np.empty(len(reader.integ)*globalMesh.NGauss*(reader.time.shape[0]))
      nGpE = globalMesh.NGaussperEl
      count0 = 0
      for timeStep in range(reader.time.shape[0]):
        field = np.empty((len(reader.integ),globalMesh.NGauss))
        count = 0
        for sd in range(self.nbsd):
          for k in range(len(reader.integ)):
            count = 0
            for el in localIdstotreat[sd]:
              field[k,el*nGpE:(el+1)*nGpE] = localDataInteg[sd][reader.integ[k]][count:count+nGpE,timeStep]
              count += nGpE
        for m in range(len(globalIdstotreat)):
            for k in range(len(reader.integ)):
              data_integ[count0:count0+nGpE] = field[k,nGpE*m:nGpE*(m+1)]
              count0 += nGpE

      data_integ.astype(np.float32).byteswap().tofile(self.outputFolder + self.name + ".integ")

      #write .node
      data_node = np.zeros(len(reader.node)*globalMesh.Nodes*(reader.time.shape[0]))
      for timeStep in range(reader.time.shape[0]):
        count = 0
        for sd in range(self.nbsd):
          for k in range(len(reader.node)):
            indices = list(map(lambda x: x + len(reader.node)*globalMesh.Nodes*timeStep+k*globalMesh.Nodes, localoriginalIDNodes[sd]))
            data_node[indices] = localDataNode[sd][reader.node[k]][:,timeStep]

      data_node.astype(np.float32).byteswap().tofile(self.outputFolder + self.name + ".node")

      #write .ut
      __string = "**meshfile " + os.path.relpath(self.dataFolder, self.outputFolder) + os.sep + cutGeof+"\n"
      with open(self.dataFolder + self.name + "-001.ut", 'r') as inFile:
        inFile.readline()
        for i in range(3):
          __string += inFile.readline()

      with open(self.outputFolder + self.name + ".ut", "w") as outFile:
        outFile.write(__string)
        for timeStep in range(reader.time.shape[0]):
          line = ""
          for i in range(4):
            line += str(int(reader.time[timeStep,i]))+" "
          line += str(reader.time[timeStep,4])+"\n"
          outFile.write(line)

def GeofFromCut(dataFolder, cutName):
  cutFile = open(dataFolder+cutName+".cut", 'r')
  strings = cutFile.readlines()
  return strings[1].split()[1]

def Tag3D(mesh):
  for name, data in mesh.elements.items():
    if EN.dimension[name] == 3:
      mesh.GetElementsOfType(name).tags.CreateTag('3D').SetIds(mesh.GetElementsOfType(name).originalIds-1)


def Return3DElements(mesh):
  for name, data in mesh.elements.items():
    if '3D' in data.tags:
      idstotreat = data.tags['3D'].GetIds()
      mesh.NnodeperEl = EN.numberOfNodes[name]
      mesh.p, mesh.w =  Lagrange(name)
      mesh.NGaussperEl = len(mesh.w)
      mesh.NGauss = data.GetNumberOfElements()*mesh.NGaussperEl
      mesh.nbElements = data.GetNumberOfElements()
      mesh.Nodes = mesh.GetNumberOfNodes()
  return idstotreat


def CheckIntegrity():

    from BasicTools.Helpers.Tests import TestTempDir
    tempdir = TestTempDir.GetTempPath()
    import BasicTools.TestData as BasicToolsTestData

    ##################################
    # EXEMPLE SYNTAXE DU MERGER
    import BasicTools.IO.Parallel.UtMerger as UM
    merger = UM.UtMerger()
    merger.SetName("cube")
    merger.SetdataFolder(BasicToolsTestData.GetTestDataPath() + "UtParExample/")
    merger.SetOutputFolder(tempdir)
    merger.Merge()
    ##################################

    import filecmp
    print(TFormat.InRed("node files equals  ? "+ str(filecmp.cmp(tempdir + "cube.node",  BasicToolsTestData.GetTestDataPath() + "UtParExample/cube.node", shallow=False))))
    print(TFormat.InRed("integ files equals ? "+ str(filecmp.cmp(tempdir + "cube.integ", BasicToolsTestData.GetTestDataPath() + "UtParExample/cube.integ", shallow=False))))
    print(tempdir)
    return "ok"

if __name__ == '__main__':
    print((CheckIntegrity()))# pragma: no cover
