# -*- coding: utf-8 -*-
""" Geof file reader (Zset mesh file)

"""
import numpy as np

__author__ = "Felipe Bordeu"

import BasicTools.Containers.ElementNames as EN

from BasicTools.IO.ReaderBase import ReaderBase

GeofNumber = {}
PermutationZSetToBasicTools = {}

GeofNumber['l2d1']   = EN.Point_1
GeofNumber['l2d2']   = EN.Bar_2
GeofNumber['quad']   = EN.Bar_3

GeofNumber['q4'] = EN.Quadrangle_4
GeofNumber['quad4'] = EN.Quadrangle_4
GeofNumber['c2d4'] = EN.Quadrangle_4
GeofNumber['c2d8'] = EN.Quadrangle_8
GeofNumber["c2d3"] = EN.Triangle_3
GeofNumber["c2d6"] = EN.Triangle_6
PermutationZSetToBasicTools["c2d6"] = [0, 2, 4, 1, 3, 5]
PermutationZSetToBasicTools["c3d6"] = [0, 2, 4, 1, 3, 5]


GeofNumber['c3d4'] = EN.Tetrahedron_4
GeofNumber['c3d4r'] = EN.Tetrahedron_4
GeofNumber['c3d6'] = EN.Wedge_6
GeofNumber['c3d8'] = EN.Hexaedron_8

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


nbIntegrationsPoints = {}
nbIntegrationsPoints["b2d3"] = 3;
nbIntegrationsPoints["b2d3r"] = -1;
nbIntegrationsPoints["bn2d3"] = 3;
nbIntegrationsPoints["bn3d4"] = 4;
nbIntegrationsPoints["bn3d6"] = 6;
nbIntegrationsPoints["bn3d8"] = 8;
nbIntegrationsPoints["c1d2"] = 3;
nbIntegrationsPoints["c1d2r"] = 1;
nbIntegrationsPoints["c1d3"] = -1;
nbIntegrationsPoints["c1d3r"] = -1;
nbIntegrationsPoints["c2d3"] = 3;
nbIntegrationsPoints["c2d3b"] = 3;
nbIntegrationsPoints["c2d3gp16"] = 16;
nbIntegrationsPoints["c2d3gp25"] = 25;
nbIntegrationsPoints["c2d3gp36"] = 36;
nbIntegrationsPoints["c2d3gp49"] = 49;
nbIntegrationsPoints["c2d3gp64"] = 64;
nbIntegrationsPoints["c2d3gp81"] = 81;
nbIntegrationsPoints["c2d3gp9"] = 9;
nbIntegrationsPoints["c2d3p3"] = 16;
nbIntegrationsPoints["c2d3r"] = 1;
nbIntegrationsPoints["c2d4"] = 4;
nbIntegrationsPoints["c2d4gp12"] = 12;
nbIntegrationsPoints["c2d4gp13"] = 13;
nbIntegrationsPoints["c2d4gp64"] = 64;
nbIntegrationsPoints["c2d4r"] = 1;
nbIntegrationsPoints["c2d6"] = 6;
nbIntegrationsPoints["c2d6gp16"] = 16;
nbIntegrationsPoints["c2d6gp25"] = 25;
nbIntegrationsPoints["c2d6gp36"] = 36;
nbIntegrationsPoints["c2d6gp4"] = 4;
nbIntegrationsPoints["c2d6gp64"] = 64;
nbIntegrationsPoints["c2d6gp81"] = 81;
nbIntegrationsPoints["c2d6gp9"] = 9;
nbIntegrationsPoints["c2d6r"] = 3;
nbIntegrationsPoints["c2d8"] = 9;
nbIntegrationsPoints["c2d8gp12"] = 12;
nbIntegrationsPoints["c2d8gp13"] = 13;
nbIntegrationsPoints["c2d8gp64"] = 64;
nbIntegrationsPoints["c2d8r"] = 4;
nbIntegrationsPoints["c2dgp49"] = -1;
nbIntegrationsPoints["c3d10"] = 5;
nbIntegrationsPoints["c3d10_11"] = 11;
nbIntegrationsPoints["c3d10_4"] = 4;
nbIntegrationsPoints["c3d10r"] = 1;
nbIntegrationsPoints["c3d12"] = 12;
nbIntegrationsPoints["c3d12f"] = 18;
nbIntegrationsPoints["c3d12l"] = -1;
nbIntegrationsPoints["c3d12r"] = 6;
nbIntegrationsPoints["c3d13"] = 27;
nbIntegrationsPoints["c3d13faux"] = 5;
nbIntegrationsPoints["c3d13r"] = 5;
nbIntegrationsPoints["c3d15"] = 18;
nbIntegrationsPoints["c3d15_9"] = 9;
nbIntegrationsPoints["c3d15l"] = 0;
nbIntegrationsPoints["c3d15r"] = 6;
nbIntegrationsPoints["c3d16"] = 18;
nbIntegrationsPoints["c3d16f"] = 27;
nbIntegrationsPoints["c3d16l"] = 0;
nbIntegrationsPoints["c3d16r"] = 8;
nbIntegrationsPoints["c3d20"] = 27;
nbIntegrationsPoints["c3d20l"] = 0;
nbIntegrationsPoints["c3d20r"] = 8;
nbIntegrationsPoints["c3d4"] = 4;
nbIntegrationsPoints["c3d4_11"] = 1;
nbIntegrationsPoints["c3d4_35"] = 1;
nbIntegrationsPoints["c3d4b"] = 4;
nbIntegrationsPoints["c3d4r"] = 1;
nbIntegrationsPoints["c3d5"] = 5;
nbIntegrationsPoints["c3d5_27"] = 27;
nbIntegrationsPoints["c3d5_6"] = 6;
nbIntegrationsPoints["c3d6"] = 6;
nbIntegrationsPoints["c3d6l"] = -1;
nbIntegrationsPoints["c3d6r"] = 2;
nbIntegrationsPoints["c3d8"] = 8;
nbIntegrationsPoints["c3d8l"] = -1;
nbIntegrationsPoints["c3d8r"] = 1;
nbIntegrationsPoints["cax3"] = 3;
nbIntegrationsPoints["cax3b"] = 3;
nbIntegrationsPoints["cax3gp16"] = 16;
nbIntegrationsPoints["cax3gp25"] = 25;
nbIntegrationsPoints["cax3gp36"] = 36;
nbIntegrationsPoints["cax3gp4"] = 4;
nbIntegrationsPoints["cax3gp49"] = 49;
nbIntegrationsPoints["cax3gp64"] = 64;
nbIntegrationsPoints["cax3gp81"] = 81;
nbIntegrationsPoints["cax3gp9"] = 9;
nbIntegrationsPoints["cax3r"] = 1;
nbIntegrationsPoints["cax4"] = 4;
nbIntegrationsPoints["cax4gp12"] = 12;
nbIntegrationsPoints["cax4gp13"] = 13;
nbIntegrationsPoints["cax4gp64"] = 64;
nbIntegrationsPoints["cax4r"] = 1;
nbIntegrationsPoints["cax6"] = 6;
nbIntegrationsPoints["cax6gp16"] = -1;
nbIntegrationsPoints["cax6gp25"] = -1;
nbIntegrationsPoints["cax6gp36"] = -1;
nbIntegrationsPoints["cax6gp49"] = -1;
nbIntegrationsPoints["cax6gp64"] = -1;
nbIntegrationsPoints["cax6gp81"] = -1;
nbIntegrationsPoints["cax6gp9"] = -1;
nbIntegrationsPoints["cax6r"] = 3;
nbIntegrationsPoints["cax8"] = 9;
nbIntegrationsPoints["cax8gp12"] = 12;
nbIntegrationsPoints["cax8gp13"] = 13;
nbIntegrationsPoints["cax8gp64"] = 64;
nbIntegrationsPoints["cax8r"] = 4;
nbIntegrationsPoints["ct2d3"] = 0;
nbIntegrationsPoints["ct3d4"] = 0;
nbIntegrationsPoints["cyl2"] = 3;
nbIntegrationsPoints["cyl2r"] = 1;
nbIntegrationsPoints["deb2"] = -1;
nbIntegrationsPoints["dummy2"] = -1;
nbIntegrationsPoints["dummy21"] = -1;
nbIntegrationsPoints["dummy22"] = -1;
nbIntegrationsPoints["dummy23"] = -1;
nbIntegrationsPoints["dummy24"] = -1;
nbIntegrationsPoints["dummy25"] = -1;
nbIntegrationsPoints["dummy26"] = -1;
nbIntegrationsPoints["dummy3"] = -1;
nbIntegrationsPoints["dummy31"] = -1;
nbIntegrationsPoints["dummy32"] = -1;
nbIntegrationsPoints["dummy33"] = -1;
nbIntegrationsPoints["dummy34"] = -1;
nbIntegrationsPoints["dummy35"] = -1;
nbIntegrationsPoints["dummy36"] = -1;
nbIntegrationsPoints["gen10"] = -1;
nbIntegrationsPoints["gen11"] = -1;
nbIntegrationsPoints["gen12"] = -1;
nbIntegrationsPoints["gen13"] = -1;
nbIntegrationsPoints["gen14"] = -1;
nbIntegrationsPoints["gen15"] = -1;
nbIntegrationsPoints["gen2"] = -1;
nbIntegrationsPoints["gen3"] = -1;
nbIntegrationsPoints["gen3dX"] = -1;
nbIntegrationsPoints["gen4"] = -1;
nbIntegrationsPoints["gen5"] = -1;
nbIntegrationsPoints["gen6"] = -1;
nbIntegrationsPoints["gen7"] = -1;
nbIntegrationsPoints["gen8"] = -1;
nbIntegrationsPoints["gen9"] = -1;
nbIntegrationsPoints["i2d4"] = 2;
nbIntegrationsPoints["i2d4r"] = 1;
nbIntegrationsPoints["i3d12"] = 6;
nbIntegrationsPoints["i3d12r"] = 3;
nbIntegrationsPoints["i3d6"] = 3;
nbIntegrationsPoints["i3d8"] = 4;
nbIntegrationsPoints["i3d8r"] = 1;
nbIntegrationsPoints["iax4"] = 2;
nbIntegrationsPoints["iax4r"] = 1;
nbIntegrationsPoints["l2d2"] = -1;
nbIntegrationsPoints["l2d2r"] = -1;
nbIntegrationsPoints["l2d3"] = 1;
nbIntegrationsPoints["l2d3r"] = 1;
nbIntegrationsPoints["l3d1"] = -1;
nbIntegrationsPoints["l3d1r"] = -1;
nbIntegrationsPoints["l3d2"] = 1;
nbIntegrationsPoints["l3d2r"] = 1;
nbIntegrationsPoints["l3d3"] = 1;
nbIntegrationsPoints["l3d3r"] = 1;
nbIntegrationsPoints["lag1"] = -1;
nbIntegrationsPoints["m3d3"] = -1;
nbIntegrationsPoints["m3d3r"] = -1;
nbIntegrationsPoints["m3d4"] = -1;
nbIntegrationsPoints["m3d4r"] = -1;
nbIntegrationsPoints["m3d6"] = -1;
nbIntegrationsPoints["m3d6r"] = -1;
nbIntegrationsPoints["m3d8"] = -1;
nbIntegrationsPoints["m3d8r"] = -1;
nbIntegrationsPoints["r2d2"] = 0;
nbIntegrationsPoints["r3d2"] = 0;
nbIntegrationsPoints["rve1d"] = -1;
nbIntegrationsPoints["rve2d"] = -1;
nbIntegrationsPoints["rve3d"] = -1;
nbIntegrationsPoints["s3d4"] = 8;
nbIntegrationsPoints["s3d6r"] = 6;
nbIntegrationsPoints["s3d8r"] = 8;
nbIntegrationsPoints["sax2"] = 3;
nbIntegrationsPoints["sax2r"] = 1;
nbIntegrationsPoints["sax2rr"] = -1;
nbIntegrationsPoints["sax3"] = -1;
nbIntegrationsPoints["sax3r"] = -1;
nbIntegrationsPoints["sax3rr"] = -1;
nbIntegrationsPoints["skin2d3"] = 3;
nbIntegrationsPoints["skin3d4"] = 4;
nbIntegrationsPoints["skin3d4r"] = 1;
nbIntegrationsPoints["skin3d6"] = 6;
nbIntegrationsPoints["skin3d6r"] = 3;
nbIntegrationsPoints["skin3d8r"] = 4;
nbIntegrationsPoints["skinax3"] = 3;
nbIntegrationsPoints["sph2"] = 3;
nbIntegrationsPoints["sph2r"] = 1;
nbIntegrationsPoints["spr1"] = 0;
nbIntegrationsPoints["spr1_2d"] = 0;
nbIntegrationsPoints["spr1_3d"] = 0;
nbIntegrationsPoints["spr2"] = -1;
nbIntegrationsPoints["x2d3"] = 3;
nbIntegrationsPoints["x2d3r"] = 1;
nbIntegrationsPoints["x2d4"] = 4;
nbIntegrationsPoints["x2d4r"] = 1;
nbIntegrationsPoints["x2d6"] = 6;
nbIntegrationsPoints["x2d6r"] = 3;
nbIntegrationsPoints["x2d8"] = 9;
nbIntegrationsPoints["x2d8r"] = 4;
nbIntegrationsPoints["x3d10"] = 5;
nbIntegrationsPoints["x3d10r"] = 1;
nbIntegrationsPoints["x3d15"] = 18;
nbIntegrationsPoints["x3d15r"] = 6;
nbIntegrationsPoints["x3d20"] = 27;
nbIntegrationsPoints["x3d20r"] = 8;
nbIntegrationsPoints["x3d4"] = 4;
nbIntegrationsPoints["x3d4r"] = 1;
nbIntegrationsPoints["x3d6"] = 6;
nbIntegrationsPoints["x3d6r"] = 2;
nbIntegrationsPoints["x3d8"] = 8;
nbIntegrationsPoints["x3d8r"] = 1;
nbIntegrationsPoints["xax3"] = 3;
nbIntegrationsPoints["xax3r"] = 1;
nbIntegrationsPoints["xax4"] = 4;
nbIntegrationsPoints["xax4r"] = 1;
nbIntegrationsPoints["xax6"] = 6;
nbIntegrationsPoints["xax6r"] = 3;
nbIntegrationsPoints["xax8"] = 9;
nbIntegrationsPoints["xax8r"] = 4;

def ReadGeof(fileName=None,string=None,out=None,readElset=True,readFaset=True,printNotRead=True):
    reader = GeofReader()
    reader.SetFileName(fileName)
    reader.SetStringToRead(string)
    reader.Read(fileName=fileName, string=string,out=out,readElset=readElset,readFaset=readFaset,printNotRead=printNotRead)
    return reader.output

class GeofReader(ReaderBase):
  def __init__(self):
        super(GeofReader,self).__init__()
        self.commentChar= "%"
        self.readFormat = 'r'

  def ReadMetaData(self):
    res = {}
    self.StartReading()
    nbIP  = 0
    l = self.ReadCleanLine()

    while(True):
      if l == None :
          break
      if l[0] != "*":
          pass
      elif l.find("**node")>-1:
        l       = self.ReadCleanLine()
        s       = l.split()
        res["nbNodes"] = int(s[0])
      elif l.find("**element")>-1:
        l  = self.ReadCleanLine()
        res['nbElements'] = int(l.split()[0])
        IPPerElement = np.empty(res['nbElements'],dtype= np.int)
        cpt = 0
        while(True):
          l  = self.ReadCleanLine()
          if l.find("**") > -1:
               break
          s = l.split()
          IPPerElement[cpt] = nbIntegrationsPoints[s[1]]
          cpt +=1
        continue
      l = self.ReadCleanLine()

    res['nbIntegrationPoints'] = np.sum(IPPerElement)
    res['IPPerElement'] = IPPerElement
    self.EndReading()

    return res


  def Read(self, fileName=None,string=None,out=None,readElset=True,readFaset=True,printNotRead=True):
    import BasicTools.Containers.UnstructuredMesh as UM

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
        res.originalIDNodes= np.empty((nbNodes,),dtype=np.int)
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
        self.PrintVerbose( "nbElements {}".format(nbElements) )
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
          if l.find("**") > -1 or readElset == False:
               break
          s = l.split()

          for soid in s:
              oid = int(soid)
              #res.AddElementToTagUsingOriginalId(int(oid),elsetname)
              oidToElementContainer[oid].tags.CreateTag(elsetname,False).AddToTag(oidToLocalElementNumber[oid])
        continue

      if l.find("**faset")>-1 or l.find("**liset")>-1:
        fasetName = l[8:]
        print("Reading Group " + fasetName)
        while(True):
          l  = self.ReadCleanLine()
          if l.find("**") > -1 or readFaset == False:
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
      if printNotRead == True:
        print("line starting with <<"+l[:20]+">> not considered in the reader")
      l = self.ReadCleanLine()
      continue

    self.EndReading()
    res.PrepareForOutput()
    self.output = res
    return res


from BasicTools.IO.IOFactory import RegisterReaderClass
RegisterReaderClass(".geof",GeofReader)


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


    import BasicTools.Containers.UnstructuredMesh as UM
    res = ReadGeof(fileName=newFileName,out = UM.UnstructuredMesh())
    print(res)
    return 'ok'

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
