# -*- coding: utf-8 -*-
""" Geof file reader (Zset mesh file)

"""
import numpy as np

__author__ = "Felipe Bordeu"

import BasicTools.FE.ElementNames as EN
from BasicTools.IO.ReaderBase import ReaderBase

GeofNumber = {}
PermutationZSetToBasicTools = {}

GeofNumber['l2d1']   = EN.Point_1
GeofNumber['l2d2']   = EN.Bar_2
GeofNumber['quad']   = EN.Bar_3

GeofNumber['quad4'] = EN.Quadrangle_4
GeofNumber['c2d4'] = EN.Quadrangle_4
GeofNumber["c2d3"] = EN.Triangle_3
GeofNumber["c2d6"] = EN.Triangle_6
PermutationZSetToBasicTools["c2d6"] = [0, 2, 4, 1, 3, 5]
PermutationZSetToBasicTools["c3d6"] = [0, 2, 4, 1, 3, 5]


GeofNumber['c3d4'] = EN.Tetrahedron_4
GeofNumber['c3d10'] = EN.Tetrahedron_10
PermutationZSetToBasicTools["c3d10"] = [2, 0, 1, 9, 5, 3, 4, 8, 6, 7]
GeofNumber['c3d20'] = EN.Hexaedron_20
PermutationZSetToBasicTools["c3d20"] = [0,2,4,6,12,14,16,18,1,3,5,7,13,15,17,19,8,9,10,11]

GeofNumber['t3']   = EN.Triangle_3

GeofNumber['t6']   = EN.Triangle_6
PermutationZSetToBasicTools["t6"] = PermutationZSetToBasicTools["c3d6"]
GeofNumber['q8']   = EN.Quadrangle_8
PermutationZSetToBasicTools["q8"] = [0, 2, 4, 6, 1, 3, 5, 7]

PermutationZSetToBasicTools["c3d15"] = [0, 2, 4, 9, 11, 13, 1, 3, 5, 10, 12, 14, 6, 7, 8]


def ReadGeof(fileName=None,string=None,out=None):
    reader = GeofReader()
    reader.SetFileName(fileName)
    reader.SetStringToRead(string)
    reader.Read(fileName=fileName, string=string,out=out)
    return reader.output

class GeofReader(ReaderBase):
  def __init__(self):
        super(GeofReader,self).__init__()
        self.commentChar= "%"
        self.readFormat = 'r'


  def Read(self, fileName=None,string=None,out=None):
    import BasicTools.FE.UnstructuredMesh as UM

    if fileName is not None:
      self.SetFileName(fileName)

    if string is not None:
      self.SetStringToRead(string)

    self.StartReading()

    if out is None:
        res = UM.UnstructuredMesh()
    else:
        res = out


    filetointernalid = {}

    oidToElementContainer = {}
    oidToLocalElementNumber = {}
    l = self.ReadCleanLine()
    while(True):

      #premature EOF
      if l is None:
          print("ERROR premature EOF: please check the integrity of your geof file") # pragma: no cover
          break # pragma: no cover


      #if len(l) == 0: l = string.readline().strip('\n').lstrip().rstrip(); continue

      if l.find("**node")>-1:
        l       = self.ReadCleanLine()
        s       = l.split()
        nbNodes = int(s[0])
        dim     = int(s[1])
        print("Reading "+str(nbNodes)+ " Nodes in dimension "+str(dim))
        res.nodes = np.empty((nbNodes,dim))
        res.originalIDNodes= np.empty((nbNodes,))
        cpt = 0
        while(True):
            l  = self.ReadCleanLine()
            if l.find("**") > -1:
                break
            s = l.split()
            oid = int(s[0])
            filetointernalid[oid] = cpt
            res.originalIDNodes[cpt] = int(s[0])
            res.nodes[cpt,:] = list(map(float,s[1:]))
            cpt += 1
        continue


      if l.find("**element")>-1:
        l  = self.ReadCleanLine()
        nbElements = int(l.split()[0])
        print( "nbElements {}".format(nbElements) )
        while(True):
          l  = self.ReadCleanLine()
          if l.find("**") > -1:
               break
          s = l.split()
          nametype = GeofNumber[s[1]]
          conn = [filetointernalid[x] for x in  map(int,s[2:]) ]
          elements = res.GetElementsOfType(nametype)
          oid = int(s[0])
          if s[1] in PermutationZSetToBasicTools:
              conn =  [conn[x] for x in PermutationZSetToBasicTools[s[1]] ]
          cpt = elements.AddNewElement(conn,oid)
          oidToElementContainer[oid] = elements
          oidToLocalElementNumber[oid] = cpt-1
        continue

      if l.find("**nset")>-1:
        nsetname = l.split()[1]
        print( "nset {}".format(nsetname) )

        tag = res.GetNodalTag(nsetname)

        while(True):
          l  = self.ReadCleanLine()
          if l.find("**") > -1:
               break
          s = l.split()
          for oid in s:
              tag.AddToTag(filetointernalid[int(oid)])
        continue

      if l.find("**elset")>-1:
        elsetname = l.split()[1]
        print( "elset {}".format(elsetname) )

        while(True):
          l  = self.ReadCleanLine()
          if l.find("**") > -1:
               break
          s = l.split()

          for soid in s:
              oid = int(soid)
              #res.AddElementToTagUsingOriginalId(int(oid),elsetname)
              oidToElementContainer[oid].tags.CreateTag(elsetname,False).AddToTag(oidToLocalElementNumber[oid])
        continue

      if l.find("**faset")>-1:
        fasetName = l[8:]
        print("Reading Group " + fasetName)
        while(True):
          l  = self.ReadCleanLine()
          if l.find("**") > -1:
               break
          s = l.split()
          nametype = GeofNumber[s[0]]
          conn = [filetointernalid[x] for x in  map(int,s[1:])]

          if s[0] in PermutationZSetToBasicTools:
              conn =  [conn[x] for x in PermutationZSetToBasicTools[s[0]] ]

          elements = res.GetElementsOfType(nametype)
          localId = elements.AddNewElement(conn,-1)
          elements.GetTag(fasetName).AddToTag(localId-1)
        continue

      if l.find("***return")>-1:
        print("End file")
        break

      if l.find("***geometry")>-1 or l.find("***group")>-1:
        l = self.ReadCleanLine()
        continue

      #case not treated
      print("line starting with <<"+l[:20]+">> not considered in the reader")
      l = self.ReadCleanLine()
      continue

    self.EndReading()
    res.PrepareForOutput()
    self.output = res
    return res



def CheckIntegrity():
    data = u"""
    ***geometry
    **node
    5 3
    1 0.0000000000000000e+00          0.0000000000000000e+00          0.0000000000000000e+00
    2 6.0000000000000019e-02          0.0000000000000000e+00          0.0000000000000000e+00
    3 6.0000000000000012e-02          7.4999999999999928e-02          0.0000000000000000e+00
    4 2.0000000000000000e-02          1.7999999999999900e-02          3.0000000000000000e-02
    5 3.0000000000000019e-02          0.0000000000000000e+00          0.0200000000000000e+00
    **element
    1
    1 c3d4 1 2 3 4
    ***group
    **nset g1
    1 2 3
    **elset g2
    1
    **faset tri
    t3  1 2 3
    **faset quads
    quad 1 2 5

    **not treated
    ***return
    """

    res = ReadGeof(string=data)
    print(res)

    from BasicTools.Helpers.Tests import WriteTempFile
    newFileName = WriteTempFile(filename="GeofFileTest.geof",content=data)


    import BasicTools.FE.UnstructuredMesh as UM
    res = ReadGeof(fileName=newFileName,out = UM.UnstructuredMesh())
    print(res)
    return 'ok'

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
