# -*- coding: utf-8 -*-
__author__ = "Felipe Bordeu"

import numpy as np

from BasicTools.Helpers.BaseOutputObject import BaseOutputObject

class Loading(object):
  def __init__(self):
    self.cycle_number = None
    self.sequence_number = None
    self.increment = None
    self.time = None
    self.temperature = {}
    self.pressure = {}

  def __str__(self):
    res  = "\n    cycle numbers       : " + str(self.cycle_number)
    res += "\n    sequence numbers    : " + str(self.sequence_number)
    res += "\n    increments sequence : " + str(self.increment)
    res += "\n    time sequence       : " + str(self.time)
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
    self.data_node  = {}
    self.data_ctnod = {}
    self.data_integ = {}


class ProblemData(BaseOutputObject):
    def __init__(self):
       super(ProblemData,self).__init__()
       self.name = ""
       self.folder = ""
       self.heading = ""
       self.orientations = {}
       self.materials = {}
       self.studyCases = {}
       self.solutions = {}
       self.loadings = {}
       self.mesh = None


    def Write(self, loadingKey):
       import BasicTools.IO.UtWriter as UW
       UtW = UW.UtWriter()
       UtW.SetName(self.name)
       UtW.SetFolder(self.folder)
       UtW.AttachMesh(self.mesh)
       UtW.AttachData(self.solutions[loadingKey].data_node, self.solutions[loadingKey].data_ctnod, self.solutions[loadingKey].data_integ)
       UtW.AttachSequence(self.loadings[loadingKey].cycle_number, self.loadings[loadingKey].sequence_number, self.loadings[loadingKey].increment, self.loadings[loadingKey].time)
       UtW.Write(writeGeof=True)


    def __str__(self):
        res = "  Name : " + self.name + "\n"
        res += "  Folder : " + self.folder + "\n"
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
        res += "\n    first : " + str(self.first)
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

    loading = Loading()
    loading.cycle_number = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]
    loading.sequence_number = [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
    loading.increment = [0, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    loading.time = [0., 0.2, 0.4, 0.6, 0.8, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.2, 2.4, 2.6, 2.8, 3., 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4., 4.2, 4.4, 4.6, 4.8, 5., 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.]
    res.loadings["Run1"] = loading

    import BasicTools.IO.UtReader as UR
    reader = UR.UtReader()
    reader.SetFileName(BasicToolsTestData.GetTestDataPath() + "UtExample/cube.ut")
    reader.ReadMetaData()

    reader.atIntegrationPoints = False
    Nnode = reader.Read(fieldname="U1", timeIndex=0).shape[0]
    reader.atIntegrationPoints = True
    Nint = reader.Read(fieldname="sig11", timeIndex=0).shape[0]
    Ntime = len(loading.time)

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

    solution = Solution()
    solution.data_node  = data_node
    solution.data_ctnod = data_ctnod
    solution.data_integ = data_integ
    res.solutions["Run1"] = solution

    import BasicTools.IO.GeofReader as GR
    mymesh = GR.ReadGeof(fileName=BasicToolsTestData.GetTestDataPath() + "UtExample/cube.geof")

    res.name   = "toto"
    res.folder = tempdir
    res.mesh   = mymesh
  
    ##################################
    # Check WriteUt  
    res.Write("Run1")
    ##################################

    print("Temp directory =", tempdir)
    import filecmp
    print("node files equals  ?", filecmp.cmp(tempdir + "toto.node",  BasicToolsTestData.GetTestDataPath() + "UtExample/toto.node", shallow=False))
    print("ctnod files equals ?", filecmp.cmp(tempdir + "toto.ctnod", BasicToolsTestData.GetTestDataPath() + "UtExample/toto.ctnod", shallow=False))
    print("integ files equals ?", filecmp.cmp(tempdir + "toto.integ", BasicToolsTestData.GetTestDataPath() + "UtExample/toto.integ", shallow=False))
    print("ut files equals    ?", filecmp.cmp(tempdir + "toto.ut",    BasicToolsTestData.GetTestDataPath() + "UtExample/toto.ut", shallow=False))


    print(res)
    return "ok"

if __name__ == '__main__':
    print((CheckIntegrity()))# pragma: no cover
