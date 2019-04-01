# -*- coding: utf-8 -*-
__author__ = "Felipe Bordeu"

import numpy as np
from BasicTools.Helpers.TextFormatHelper import TFormat

from BasicTools.Helpers.BaseOutputObject import BaseOutputObject
from BasicTools.FE.Fields.ConstantField import ConstantField
from BasicTools.FE.Fields.FEField import FEField
from BasicTools.FE.Fields.IntegrationPointField import IntegrationPointField

class Loading(object):
  def __init__(self):
    self.timeSequence = None

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
       #loading.Ntime = loading.timeSequence.shape[0]
       self.loadings[tag] = loading

    def AttachSolution(self, tag, solution):
       self.solutions[tag] = solution

    def InitWriter(self, loadingKey, name, folder, Nnode = None, Nint = None):
       import BasicTools.IO.UtWriter as UW
       import os
       self.UtW = UW.UtWriter()
       self.UtW.SetFileName(folder+os.sep+name)
       self.UtW.AttachMesh(self.mesh)
       self.UtW.AttachDataFromProblemData(self, loadingKey, Nnode = Nnode, Nint = Nint)
       self.UtW.AttachSequence(self.loadings[loadingKey].timeSequence)

    def Write(self, loadingKey, name, folder, skipCtnod = False):
       self.InitWriter(loadingKey, name, folder)
       self.UtW.WriteFiles(writeGeof=True, skipCtnod = skipCtnod)


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


class Transform(BaseOutputObject):
    def __init__(self, offset=None, first=None, second=None):
       super(Transform,self).__init__()
       self.system = ''

       self.offset  = np.array([0.0, 0.0, 0.0], dtype=np.float)
       self.RMatrix = np.array([[1.0,0,0],[0,1,0],[0,0,1]])
       self.keepOrthogonal = True
       self.keepNormalised = True

       # offset off the new origin with respect to the old
       if offset is None:
           self.SetOffset([0.0, 0.0, 0.0])
       else:
           self.SetOffset(offset)

       if first is None:
           self.SetFirst([1.0, 0.0, 0.0])
       else:
           self.SetFirst(first)

       if second is None:
           self.SetSecond([0.0, 1.0, 0.0])
       else:
           self.SetSecond(second)

    def GetDirection(self,i, pos=None, direction=None):
        return self.RMatrix[i,:]

    def SetOffset(self,data):
        # this point define the origin of the new coordinate system
        self.offset  = np.array(data, dtype=np.float)

    def SetFirst(self,data):
        # this point define the x coordinate (direction) with respect to the new origin
        first  = np.array(data, dtype=np.float)
        if self.keepNormalised :
            first /= np.linalg.norm(first)
        self.RMatrix[0,:] = first
        self.SetSecond(self.RMatrix[1,:])
        self.SetThird(self.RMatrix[2,:])

    def SetSecond(self,data):
        # this point define the y coordinate (direction) with respect to the new origin
        # the z direction is calculated with a cross product
        first = self.RMatrix[0,:]
        second = np.array(data, dtype=np.float)
        second -= first*np.dot(first,second)/np.dot(first,first)
        if self.keepNormalised :
            second /= np.linalg.norm(second)
        self.RMatrix[1,:] = second
        self.SetThird(self.RMatrix[2,:])

    def SetThird(self,data=None):
        first = self.RMatrix[0,:]
        second = self.RMatrix[1,:]

        if data is None:
            self.RMatrix[2,:] = np.cross(first, second)
        else:
            third = np.array(data, dtype=np.float)
            third -= first*np.dot(first,third)/np.dot(first,first)
            third -= second*np.dot(second,third)/np.dot(second,second)
            self.RMatrix[2,:] = third

        if self.keepNormalised :
            self.RMatrix[2,:] /= np.linalg.norm(self.RMatrix[2,:])

    def SetOperator(self, first=None, second=None, third=None, op=None)    :
        if op in None:
            self.SetFirst(first)
            self.SetSecond(second)
            self.SetThird(third)
        else:
            if (first is not None) or (second is not None) or  (third is not None):
                raise(Exception("Cant define operanter and direction verctor at the same time") )
            self.RMatrix[:,:] = op[:,:]

            if self.keepNormalised :
                for i in [0,1,2]:
                    self.RMatrix[i,:] /= np.linalg.norm(self.RMatrix[i,:])

    def ApplyInvTransform(self,point):
        # we apply inverse of the transformation
        #p = point+self.offset
        if self.keepNormalised == False:
            return np.dot(np.linalg.inv(self.RMatrix),point)+self.offset
        else:
            return np.dot(self.RMatrix.T,point)+self.offset

    def ApplyTransform(self,point):
        # we apply inverse of the transformation
        #p = point+self.offset
        return np.dot(self.RMatrix,point-self.offset)

    def GetOrthoNormalBase(self):
        return Transform(self.offset, self.RMatrix[0,:], self.RMatrix[1,:])

    def __str__(self):
        res  = "\n    offset : " + str(self.offset)
        res += "\n    first  : " + str(self.RMatrix[0,:])
        res += "\n    second : " + str(self.RMatrix[1,:])
        res += "\n    third : " + str(self.RMatrix[2,:])
        return res


def CheckIntegrity(GUI=False):

    res = ProblemData()

    orient = Transform()
    orient.keepNormalised = False
    orient.SetFirst([2,1,0])
    orient.SetSecond([1,3,0])
    p = [1,3,3]
    print(orient)

    pprim = orient.ApplyTransform(p)
    print("p",p)
    print("pprim",pprim)

    print("p back",orient.ApplyInvTransform(pprim))
    if np.linalg.norm(orient.ApplyInvTransform(pprim) -p ) > 1e-15 :raise

    res.orientations["Toto"] = orient

    from BasicTools.Helpers.Tests import TestTempDir
    tempdir = TestTempDir.GetTempPath()
    import BasicTools.TestData as BasicToolsTestData

    import BasicTools.IO.UtReader as UR
    reader = UR.UtReader()
    reader.SetFileName(BasicToolsTestData.GetTestDataPath() + "UtExample/cube.ut")
    reader.ReadMetaData()

    reader.atIntegrationPoints = False
    Nnode = reader.ReadField(fieldname="U1", timeIndex=0).shape[0]
    reader.atIntegrationPoints = True
    Nint = reader.ReadField(fieldname="sig11", timeIndex=0).shape[0]

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
        solution.data_node[dn][:,i] = reader.ReadField(fieldname=dn, timeIndex=i)
      for dc in solution.data_ctnod:
        solution.data_ctnod[dc][:,i] = reader.ReadField(fieldname=dc, timeIndex=i)
    reader.atIntegrationPoints = True
    for i in range(Ntime):
      for di in solution.data_integ:
        solution.data_integ[di][:,i] = reader.ReadField(fieldname=di, timeIndex=i)
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
    print((CheckIntegrity(GUI=True)))# pragma: no cover
