# -*- coding: utf-8 -*-

import numpy as np
import BasicTools.FE.ElementNames as EN

#  for ABAQUS input file

InpNumber = {}

InpNumber['C3D4'] = EN.Tetrahedron_4
InpNumber['S3'] = EN.Triangle_3





def ReadInp(fileName=None,string=None):
    from cStringIO import StringIO
    import BasicTools.FE.UnstructuredMesh as UM

    if fileName is not None:
        f = open(fileName, 'r')
        string = f.read()
        f.close()

    string = StringIO(string)

    coef = 1.
    res = UM.UnstructuredMesh()
    filetointernalid = {}
    l = string.readline().strip('\n').lstrip().rstrip()

    while(True):

      if l.find("**LENGTH UNITS")>-1:
          if l.find("mm")>-1:
              coef = 0.001
          l = string.readline().strip('\n').lstrip().rstrip()
          continue

      if l.find("*NODE")>-1:
            cpt = 0
            res.originalIDNodes = np.empty((0,1), int)
            res.nodes = np.empty((0,3), float)
            l  = string.readline().strip('\n').lstrip().rstrip()
            while(True):
                if len(l) == 0: l = string.readline().strip('\n').lstrip().rstrip(); continue
                if l.find("*") > -1 or not l:
                    break
                s = l.replace(',', '').split()
                oid = int(s[0])
                if oid not in filetointernalid:
                  filetointernalid[oid] = cpt
                  res.originalIDNodes = np.vstack((res.originalIDNodes,int(s[0])))   
                  res.nodes = np.vstack((res.nodes,map(float,s[1:])))
                  cpt += 1

                l = string.readline().strip('\n').lstrip().rstrip()
            continue


      if l.find("*ELEMENT")>-1:
        s = l.replace(',', '').split()
        nametype = InpNumber[s[1][5:]]
        elset = s[2][6:]
        l = string.readline().strip('\n').lstrip().rstrip()
        cpt = 0
        while(True):
          if l.find("*") > -1 or not l:
               break
          s = l.replace(',', '').split()
          conn = [filetointernalid[x] for x in  map(int,s[1:]) ]
          elements = res.GetElementsOfType(nametype)
          oid = int(s[0])
          elements.AddNewElement(conn,oid)
          elements.GetTag(elset).AddToTag(cpt)
          cpt += 1
          l = string.readline().strip('\n').lstrip().rstrip()
        continue


      if not l: break

      #case not treated
      print("line starting with <<"+l[:20]+">> not considered in the reader")
      l = string.readline().strip('\n').lstrip().rstrip()
      continue

    res.originalIDNodes = np.squeeze(res.originalIDNodes)

    return res



def CheckIntegrity():
    data = """** -------------------------------------------------------
    ** ABAQUS input file
    ** File exported by VISUAL-ENVIRONMENT : Version10.7
    **      on 2016-4-8, at 10Hr:44min
    ** -------------------------------------------------------
    **LENGTH UNITS: mm
    *NODE
         1,      3.0310898125000696,      1.7500003247594429,  0.00059999999999149622
         2,      2.9785080000000002,     0.87500000000000011,     0.75012199999999996
         3,      3.1961719999999998,      1.0006690000000003,      1.5926899999999999
         4,      2.2943020000000001,               0.4799735,      1.4870170000000003
         5,      2.9259249999999999,                      0.,      1.4996439999999893
         6,      2.7781790000000002,     0.78753359999999994,      2.3200409999999998
    *ELEMENT, TYPE=C3D4, ELSET=AM1_labo
         1,         1,         2,         3,         4
         2,         3,         5,         6,         4"""

    res = ReadInp(string=data)

    from BasicTools.Helpers.Tests import TestTempDir
    tempdir = TestTempDir.GetTempPath()
    f =open(tempdir+"test.inp","w")
    f.write(data)
    f.close()
    res = ReadInp(fileName=tempdir+"test.inp")





    return 'ok'

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
