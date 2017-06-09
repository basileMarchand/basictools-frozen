# -*- coding: utf-8 -*-

import numpy as np
import BasicTools.FE.ElementNames as EN

GeofNumber = {}

GeofNumber['l2d1']   = EN.Point_1
GeofNumber['l2d2']   = EN.Bar_2
GeofNumber['quad4'] = EN.Quadrangle_4
GeofNumber['c2d4'] = EN.Quadrangle_4
GeofNumber["c2d3"] = EN.Triangle_3

GeofNumber['c3d4'] = EN.Tetrahedron_4
GeofNumber['c3d20'] = EN.Hexaedron_20
GeofNumber['t3']   = EN.Triangle_3
GeofNumber['q8']   = EN.Quadrangle_8


def ReadGeof(fileName=None,string=None):
    from io import StringIO
    import BasicTools.FE.UnstructuredMesh as UM

    if fileName is not None:
        string = open(fileName, 'r')
    elif string is not None:
        from io import StringIO
        string = StringIO(string)

    res = UM.UnstructuredMesh()
    filetointernalid = {}
    l = string.readline().strip('\n').lstrip().rstrip()

    while(True):
      if len(l) == 0: l = string.readline().strip('\n').lstrip().rstrip(); continue

      if l.find("**node")>-1:
            l       = string.readline().strip('\n').lstrip().rstrip()
            s       = l.split()
            nbNodes = int(s[0])
            dim     = int(s[1])
            print("Reading "+str(nbNodes)+ " Nodes in dimension "+str(dim))
            res.nodes = np.empty((nbNodes,dim))
            res.originalIDNodes= np.empty((nbNodes,))
            cpt = 0
            l  = string.readline().strip('\n').lstrip().rstrip()
            while(True):
                if len(l) == 0: l = string.readline().strip('\n').lstrip().rstrip(); continue
                if l.find("**") > -1:
                    break
                s = l.split()
                oid = int(s[0])
                filetointernalid[oid] = cpt
                res.originalIDNodes[cpt] = int(s[0])
                res.nodes[cpt,:] = list(map(float,s[1:]))
                cpt += 1
                l = string.readline().strip('\n').lstrip().rstrip()
            continue


      if l.find("**element")>-1:
        l = string.readline().strip('\n').lstrip().rstrip()
        nbElements = int(l.split()[0])
        print( "nbElements {}".format(nbElements) )
        l = string.readline().strip('\n').lstrip().rstrip()
        while(True):
          if len(l) == 0: l = string.readline().strip('\n').lstrip().rstrip(); continue
          if l.find("**") > -1:
               break
          s = l.split()
          nametype = GeofNumber[s[1]]
          conn = [filetointernalid[x] for x in  map(int,s[2:]) ]
          elements = res.GetElementsOfType(nametype)
          oid = int(s[0])
          if nametype == EN.Quadrangle_8:
              conn =  [conn[x] for x in [0, 2, 4, 6, 1, 3, 5, 7]]
          if nametype == EN.Hexaedron_20:
              conn =  [conn[x] for x in [0,2,4,6,12,14,16,18,1,3,5,7,13,15,17,19,8,9,10,11]]
          elements.AddNewElement(conn,oid)
          l = string.readline().strip('\n').lstrip().rstrip()
        continue

      if l.find("**nset")>-1:
        nsetname = l.split()[1]
        print( "nset {}".format(nsetname) )

        tag = res.GetNodalTag(nsetname)

        l = string.readline().strip('\n').lstrip().rstrip()
        while(True):
          if len(l) == 0: l = string.readline().strip('\n').lstrip().rstrip(); continue
          if l.find("**") > -1:
               break
          s = l.split()
          for oid in s:
              tag.AddToTag(filetointernalid[int(oid)])
          l = string.readline().strip('\n').lstrip().rstrip()

        continue

      if l.find("**elset")>-1:
        elsetname = l.split()[1]
        print( "elset {}".format(elsetname) )


        l = string.readline().strip('\n').lstrip().rstrip()
        while(True):
          if len(l) == 0: l = string.readline().strip('\n').lstrip().rstrip(); continue
          if l.find("**") > -1:
               break
          s = l.split()
          for oid in s:
              res.AddElementToTagUsingOriginalId(int(oid),elsetname)
          l = string.readline().strip('\n').lstrip().rstrip()

        continue

      if l.find("**faset")>-1:
        fasetName = l[8:]
        print("Reading Group " + fasetName)
        l = string.readline().strip('\n').lstrip().rstrip()
        while(True):
          if len(l) == 0: l = string.readline().strip('\n').lstrip().rstrip(); continue
          if l.find("**") > -1:
               break
          s = l.split()
          nametype = GeofNumber[s[0]]
          conn = [filetointernalid[x] for x in  map(int,s[1:])]
          if nametype == EN.Quadrangle_8:
              conn =  [conn[x] for x in [0, 2, 4, 6, 1, 3, 5, 7]]
          elements = res.GetElementsOfType(nametype)
          localId = elements.AddNewElement(conn,-1)
          elements.GetTag(fasetName).AddToTag(localId-1)
          l = string.readline().strip('\n').lstrip().rstrip()
        continue

      if l.find("***return")>-1:
        print("End file")
        break

      if l.find("**nset")>-1 or l.find("**elset")>-1:
          print(l+" not read")
          while(True):
              l = string.readline().strip('\n').lstrip().rstrip()
              if l.find("**") > -1:
                  break
          continue

      if l.find("***geometry")>-1 or l.find("***group")>-1:
        l = string.readline().strip('\n').lstrip().rstrip()
        continue

      #case not treated
      print("line starting with <<"+l[:20]+">> not considered in the reader")
      l = string.readline().strip('\n').lstrip().rstrip()
      continue

    res.PrepareForOutput()

    return res



def CheckIntegrity():
    data = u"""
    ***geometry
    **node
    4 3
    1 0.0000000000000000e+00          0.0000000000000000e+00          0.0000000000000000e+00
    2 6.0000000000000019e-02          0.0000000000000000e+00          0.0000000000000000e+00
    3 6.0000000000000012e-02          7.4999999999999928e-02          0.0000000000000000e+00
    4 0.0000000000000000e+00          7.4999999999999900e-02          0.0000000000000000e+00
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
    **not treaged
    ***return
    """

    res = ReadGeof(string=data)

    from BasicTools.Helpers.Tests import TestTempDir
    newFileName = TestTempDir().GetTempPath()+"GeofFile"
    open(newFileName,'w').write(data)
    res = ReadGeof(fileName=newFileName)
    print(res)
    return 'ok'

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
