# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#


import numpy as np
import os

from BasicTools.IO.WriterBase import WriterBase as WriterBase
import BasicTools.IO.GeofReader as GR
import BasicTools.IO.UtReader as UR
import BasicTools.Containers.ElementNames as EN
from BasicTools.FE.IntegrationsRules import LagrangeIsoParam
from BasicTools.Helpers.TextFormatHelper import TFormat

class UtSplitter(WriterBase):
    "This class can plit .ut, .goef, .ctnod, .node, .integ files frpù a monolithic Zset solution"
    def __init__(self):
        super(UtSplitter,self).__init__()
        self.timeSteps = "all"

    def __str__(self):
        res  = 'UtSplit : \n'
        res += '   Name : '+str(self.name)+'\n'
        return res

    def SetName(self,name):
        self.name = name

    def SetdataFolder(self,folder):
        self.dataFolder = folder
        from os import listdir
        from os.path import isfile, join
        temp = [f.split(".") for f in listdir(self.dataFolder + self.name + "-pmeshes" + os.sep) if isfile(join(self.dataFolder + self.name + "-pmeshes" + os.sep, f))]
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
          if 'geof' in f and self.name in f:
            for fi in f:
              if type(fi) == int:
                subdomains.append(fi)
        self.nbsd = np.max(subdomains)
        print("Number of found subdomains:", self.nbsd)

    def SetOutputFolder(self,outputFolder):
        self.outputFolder = outputFolder

    def SetTimeSteps(self,iterator):
        self.timeSteps = iterator

    def Split(self):

      reader = UR.UtReader()
      reader.SetFileName(self.dataFolder + self.name + ".ut")
      reader.ReadMetaData()


      if self.timeSteps != "all":
        reader.time = reader.time[self.timeSteps,:]
        if len(reader.time.shape) == 1:
          reader.time.shape = (1,reader.time.shape[0])


      reader.atIntegrationPoints = True
      globaldataInteg = {}
      for din in reader.integ:
        globaldataInteg[din] = np.empty((reader.meshMetadata['nbIntegrationPoints'],reader.time.shape[0]))
        for timeStep in range(reader.time.shape[0]):
          globaldataInteg[din][:,timeStep] = reader.ReadField(fieldname=din, timeIndex=int(reader.time[timeStep,0])-1)

      reader.atIntegrationPoints = False
      globaldataNode = {}
      for din in reader.node:
          globaldataNode[din] = np.empty((reader.meshMetadata['nbNodes'],reader.time.shape[0]))
          for timeStep in range(reader.time.shape[0]):
            globaldataNode[din][:,timeStep] = reader.ReadField(fieldname=din, timeIndex=int(reader.time[timeStep,0])-1)


      globalMesh = GR.ReadGeof(fileName = self.dataFolder + self.name + ".geof",readElset=False,readFaset=False,printNotRead=False)
      Tag3D(globalMesh)
      Return3DElements(globalMesh)
      globalMesh.originalIDNodes = np.array(globalMesh.originalIDNodes-1, dtype=int)

      nGpE = globalMesh.NGaussperEl


      from BasicTools.Helpers.ProgressBar import printProgressBar
      printProgressBar(0, self.nbsd, prefix = 'Progress:', suffix = 'Complete', length = 50)

      for sd in range(1,self.nbsd+1):
        sdString = "-" + str(sd).zfill(3)
        localMesh = GR.ReadGeof(fileName = self.dataFolder+self.name + "-pmeshes" + os.sep + self.name + sdString + ".geof",readElset=False,readFaset=False,printNotRead=False)
        Tag3D(localMesh)
        localIdstotreat = Return3DElements(localMesh)
        localOriginalIDNodes = np.array(localMesh.originalIDNodes-1, dtype=int)


        #write .integ
        data_integ = np.empty(len(reader.integ)*localMesh.NGauss*(reader.time.shape[0]))

        count0 = 0
        for timeStep in range(reader.time.shape[0]):
          field = np.empty((len(reader.integ),localMesh.NGauss))
          count = 0
          for k in range(len(reader.integ)):
            count = 0
            for el in localIdstotreat:
              field[k,count:count+nGpE] = globaldataInteg[reader.integ[k]][el*nGpE:(el+1)*nGpE,timeStep]
              count += nGpE
          for m in range(len(localIdstotreat)):
            for k in range(len(reader.integ)):
              data_integ[count0:count0+nGpE] = field[k,nGpE*m:nGpE*(m+1)]
              count0 += nGpE
        data_integ.astype(np.float32).byteswap().tofile(self.outputFolder + self.name + sdString + ".integ")


        #write .node
        count0 = 0
        data_node = np.zeros(len(reader.node)*localMesh.Nodes*(reader.time.shape[0]))
        for timeStep in range(reader.time.shape[0]):
          for k in range(len(reader.node)):
            data_node[count0:count0+localMesh.Nodes] = globaldataNode[reader.node[k]][localOriginalIDNodes,timeStep]
            count0 += localMesh.Nodes
        data_node.astype(np.float32).byteswap().tofile(self.outputFolder + self.name + sdString +".node")


        #write .ut
        __string = "**meshfile " + os.path.relpath(self.dataFolder, self.outputFolder) + os.sep + self.name + "-pmeshes" + os.sep + self.name + sdString + ".geof\n"
        with open(self.dataFolder + self.name + ".ut", 'r') as inFile:
          #next(inFile)
          inFile.readline()
          for i in range(3):
            __string += inFile.readline()

        with open(self.outputFolder + self.name + sdString + ".ut", "w") as outFile:
          outFile.write(__string)
          for timeStep in range(reader.time.shape[0]):
            line = ""
            for i in range(4):
              line += str(int(reader.time[timeStep,i]))+" "
            line += str(reader.time[timeStep,4])+"\n"
            outFile.write(line)

        printProgressBar(sd, self.nbsd, prefix = 'Progress:', suffix = 'Complete', length = 50)



      #write .cut
      __string = "***decomposition\n" + "  **global_mesh " + os.path.relpath(self.dataFolder, self.outputFolder) + os.sep + self.name + ".geof\n" + "  **domains "+str(self.nbsd)+"\n"
      with open(self.outputFolder + self.name + ".cut", "w") as outFile:
        outFile.write(__string)



def Tag3D(mesh):
  for name, data in mesh.elements.items():
    if EN.dimension[name] == 3:
      mesh.GetElementsOfType(name).tags.CreateTag('3D').SetIds(mesh.GetElementsOfType(name).originalIds-1)


def Return3DElements(mesh):
  for name, data in mesh.elements.items():
    if '3D' in data.tags:
      idstotreat = data.tags['3D'].GetIds()
      mesh.NnodeperEl = EN.numberOfNodes[name]
      mesh.p, mesh.w =  LagrangeIsoParam[name]
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
    # EXEMPLE SYNTAXE DU SPLITTER
    import BasicTools.IO.Parallel.UtSplitter as US
    splitter = US.UtSplitter()
    splitter.SetName("cube")
    splitter.SetdataFolder(BasicToolsTestData.GetTestDataPath() + "UtParExample/")
    splitter.SetOutputFolder(tempdir)
    splitter.Split()
    ##################################

    print("tempdir =", tempdir)

    import filecmp
    for i in range(splitter.nbsd):
      sdString = "-" + str(i+1).zfill(3)
      print(TFormat.InRed("node files for subdomain "+str(i+1)+" equals  ? "+ str(filecmp.cmp(tempdir + "cube"+sdString+".node",  BasicToolsTestData.GetTestDataPath() + "UtParExample/cube"+sdString+".node", shallow=False))))
      print(TFormat.InRed("integ files for subdomain "+str(i+1)+" equals ? "+ str(filecmp.cmp(tempdir + "cube"+sdString+".integ", BasicToolsTestData.GetTestDataPath() + "UtParExample/cube"+sdString+".integ", shallow=False))))

    print(tempdir)
    return "ok"

if __name__ == '__main__':
    print((CheckIntegrity()))# pragma: no cover
