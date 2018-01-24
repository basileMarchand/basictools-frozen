# -*- coding: utf-8 -*-
__author__ = "Felipe Bordeu"

import numpy as np
from BasicTools.Helpers.TextFormatHelper import TFormat

from BasicTools.Helpers.BaseOutputObject import BaseOutputObject
from BasicTools.FE.Fields.ConstantField import ConstantField
from BasicTools.FE.Fields.NodalField import NodalField
from BasicTools.FE.Fields.IntegrationPointField import IntegrationPointField

class Loading(object):
  def __init__(self):
    self.timeSequence = None
    self.dataLoading = {}

  def __str__(self):
    res =  TFormat.GetIndent()+TFormat.InBlue('Loading :\n')
    TFormat.II()
    res += TFormat.GetIndent()+TFormat.InBlue("Parameters: ")
    for l in self.__dict__:
      res += TFormat.GetIndent()+TFormat.InRed(l.__str__()+"  ")
    TFormat.DI()
    return res

class BoundaryConditions(object):
  def __init__(self):
    self.on = None

class StudyCase(object):
    def __init__(self):
        self.type = None

    def __str__(self):
        res = ""
#        res += str(self.name) + "\n"
        return res

class Solution(object):
    def __init__(self):
      import collections
      self.data_node  = collections.OrderedDict()
      self.data_ctnod = collections.OrderedDict()
      self.data_integ = collections.OrderedDict()

    def __str__(self):
      res =  "\n    node  : "+str(list(self.data_node.keys()))
      res += "\n    integ : "+str(list(self.data_integ.keys()))
      return res
    

class ProblemData(BaseOutputObject):
    def __init__(self):
       super(ProblemData,self).__init__()
       self.name = ""
       self.heading = ""
       self.orientations = {}
       self.materials = {}
       self.studyCases = {}
       self.solutions = {}
       self.loadings = {}
       self.mesh = None

    def AttachMesh(self, mesh):
       self.mesh = mesh

    def AttachLoading(self, tag, loading):
       loading.Ntime = loading.timeSequence.shape[0]
       self.loadings[tag] = loading

    def AttachSolution(self, tag, solution):
       self.solutions[tag] = solution

    def Write(self, loadingKey, name, folder):
       import BasicTools.IO.UtWriter as UW
       UtW = UW.UtWriter()
       UtW.SetName(name)
       UtW.SetFolder(folder)
       UtW.AttachMesh(self.mesh)
       UtW.AttachDataFromProblemData(self, loadingKey)
       UtW.AttachSequence(self.loadings[loadingKey].timeSequence)
       UtW.Write(writeGeof=True)


    def __str__(self):
        res =  "  Name : " + self.name + "\n"
        res += "  Heading : " + self.heading + "\n"
        res += "  Mesh : " + str(self.mesh)
        for l in self.loadings:
            res += "\n  Loading  '" + l + "' :"
            res +=  self.loadings[l].__str__()
        for o in self.orientations:
            res += "\n  Orientation  '" + o + "' :"
            res +=  self.orientations[o].__str__()
        for m in self.materials:
            res += "\n  Material  '" + m + "' :"
            res +=  self.materials[m].__str__()
        for s in self.solutions:
            res += "\n  Solution  '" + s + "' :"
            res +=  self.solutions[s].__str__()
        return res


class Property(object):
    def __init__(self):
        self.type = None
        self.subtype = None
        self.params = []
    def __str__(self):
        res = ""
        res += str(self.type) + " :: " + str(self.subtype) + "\n"
        res += str(self.params)
        return res


class Material(BaseOutputObject):

    #type: ELASTIC
    #subtype : ISOTROPIC
    #params (12000,0.3)

    #type: DENSITY
    #subtype : NONE
    #params (0.1,)

    def __init__(self):
        super(Material,self).__init__()
        self.name = "None"
        self.props = []

    def AddProperty(self,intype,insubtype,params):
        data = Property()
        data.type = intype
        data.subtype= insubtype
        data.params = params
        self.props.append(data)

    def __str__(self):
        res = ""
        for p in self.props:
            res += str(p) +  "\n"
        return res


class Section():
    def __init__(self):
        self.material= None
        self.elemtag = None


class Orientation(BaseOutputObject):
    def __init__(self):
       super(Orientation,self).__init__()
       self.system = ''
       # offset off the new origin with respect to the old
       self.offset = np.array([0.0, 0.0, 0.0], dtype=np.float)
       # first lies on the x axis for cartesian
       self.first  = np.array([1.0, 0.0, 0.0], dtype=np.float)
       # second lies on the y axis for cartesian
       self.second = np.array([0.0, 1.0, 0.0], dtype=np.float)

    def SetOffset(self,data):
        self.offset  = np.array(data, dtype=np.float)

    def SetFirst(self,data):
        self.first  = np.array(data, dtype=np.float)
        self.first /= np.linalg.norm(self.first)

    def SetSecond(self,data):
        self.second = np.array(data, dtype=np.float)
        self.second -= self.first*np.dot(self.first,self.second)
        self.second /= np.linalg.norm(self.second)
    def __str__(self):
        res  = "\n    offset : " + str(self.offset)
        res += "\n    first  : " + str(self.first)
        res += "\n    second : " + str(self.second)
        return res


def CheckIntegrity():
    res = ProblemData()

    orient = Orientation()
    orient.SetFirst([2,1,0])
    orient.SetSecond([1,2,0])
    res.orientations["Toto"] = orient


    from BasicTools.Helpers.Tests import TestTempDir
    tempdir = TestTempDir.GetTempPath()
    import BasicTools.TestData as BasicToolsTestData

    import BasicTools.IO.UtReader as UR
    reader = UR.UtReader()
    reader.SetFileName(BasicToolsTestData.GetTestDataPath() + "UtExample/cube.ut")
    reader.ReadMetaData()

    reader.atIntegrationPoints = False
    Nnode = reader.Read(fieldname="U1", timeIndex=0).shape[0]
    reader.atIntegrationPoints = True
    Nint = reader.Read(fieldname="sig11", timeIndex=0).shape[0]

    Ntime = reader.time.shape[0]
    NnodeVar = len(reader.node)
    NintVar = len(reader.integ)

    loading = Loading()
    loading.timeSequence = reader.time
    res.AttachLoading("Run1", loading)

    solution = Solution()
    for dnn in reader.node:
      solution.data_node[dnn] = np.empty((Nnode,Ntime))
    for din in reader.integ:
      solution.data_ctnod[din] = np.empty((Nnode,Ntime))
      solution.data_integ[din] = np.empty((Nint,Ntime))
    reader.atIntegrationPoints = False
    for i in range(Ntime):
      for dn in solution.data_node:
        solution.data_node[dn][:,i] = reader.Read(fieldname=dn, timeIndex=i)
      for dc in solution.data_ctnod:
        solution.data_ctnod[dc][:,i] = reader.Read(fieldname=dc, timeIndex=i)
    reader.atIntegrationPoints = True
    for i in range(Ntime):
      for di in solution.data_integ:
        solution.data_integ[di][:,i] = reader.Read(fieldname=di, timeIndex=i)
    res.AttachSolution("Run1", solution)

    import BasicTools.IO.GeofReader as GR
    mymesh = GR.ReadGeof(fileName=BasicToolsTestData.GetTestDataPath() + "UtExample/cube.geof")

    res.AttachMesh(mymesh)
  
    ##################################
    # Check WriteUt  
    res.Write("Run1", "toto", tempdir)
    ##################################

    print("Temp directory =", tempdir)
    import filecmp
    print("node files equals  ?", filecmp.cmp(tempdir + "toto.node",  BasicToolsTestData.GetTestDataPath() + "UtExample/cube.node", shallow=False))
    print("ctnod files equals ?", filecmp.cmp(tempdir + "toto.ctnod", BasicToolsTestData.GetTestDataPath() + "UtExample/cube.ctnod", shallow=False))
    print("integ files equals ?", filecmp.cmp(tempdir + "toto.integ", BasicToolsTestData.GetTestDataPath() + "UtExample/cube.integ", shallow=False))


    print(res)
    return "ok"

if __name__ == '__main__':
    print((CheckIntegrity()))# pragma: no cover
