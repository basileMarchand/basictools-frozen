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
    self.on = None


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

    def WriteUt(self):
       import BasicTools.IO.UtWriter as UW
       UtW = UW.UtWriter()
       for l in self.loadings:
         UtW.SetName(self.name)
         UtW.SetFolder(self.folder)
         UtW.AttachMesh(self.mesh)
         UtW.AttachData(self.solutions[l].data_node, self.solutions[l].data_ctnod, self.solutions[l].data_integ, self.solutions[l].data_node_names, self.solutions[l].data_integ_names)
         UtW.AttachSequence(cycle_number, sequence_number, increment, time)
         UtW.Write(writeGeof=True)

    def __str__(self):
        res = ""
        res = "Name : " + self.name + "\n"
        res = "Folder : " + self.folder + "\n"
        res = "Heading : " + self.heading + "\n"
        res = "Mesh : " + str(self.mesh) + "\n"
        for l in self.loadings:
            res += "Loading  '" + l + "' :\n"
            res +=  self.loadings[l].__str__()
        for o in self.orientations:
            res += "Orientation  '" + o + "' :\n"
            res +=  self.orientations[o].__str__()
        for m in self.materials:
            res += "Material  '" + m + "' :\n"
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
        res = ''
        res += "offset : " + str(self.offset) + "\n"
        res += "first : " + str(self.first) + "\n"
        res += "second : " + str(self.second) + "\n"
        return res


def CheckIntegrity():
    res = ProblemData()

    orient = Orientation()
    orient.SetFirst([2,1,0])
    orient.SetSecond([1,2,0])
    res.orientations["Toto"] = orient
    print(res)
    return "ok"

if __name__ == '__main__':
    print((CheckIntegrity()))# pragma: no cover
