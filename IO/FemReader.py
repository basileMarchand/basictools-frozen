# -*- coding: utf-8 -*-
""" Fem file reader

"""

import numpy as np

__author__ = "Felipe Bordeu"
import BasicTools.FE.UnstructuredMesh as UM
from BasicTools.IO.ReaderBase import ReaderBase
import BasicTools.FE.ElementNames as ElementNames


def ReadFem(fileName=None,string=None):
    obj = FemReader()
    if fileName is not None:
        obj.SetFileName(fileName)

    if string is not None:
        obj.SetStringToRead(string)

    return obj.Read()[0]


class FemReader(ReaderBase):
    def __init__(self,fileName = None):
        super(FemReader,self).__init__()
        self.commentChar = "$"

    def Read(self):
        self.StartReading()

        res = UM.UnstructuredMesh()
        resdic = {}
        ids = []
        xs =[]
        ys = []
        zs = []
        filetointernalid = {}
        cpt =0

        while(True):
            line = self.ReadCleanLine()
            if line is None:
                break

            l = line.split()
            key = l[0]

            if key == 'DESGLB' : continue
            if key == 'SUBCASE' : continue
            if key == 'ANALYSIS' : continue
            if key == 'SPC' : continue
            if key == 'LOAD' : continue
            if line[0:6] == 'DESOBJ' : continue
            if key == 'LABEL' : continue
            if line == "BEGIN BULK" :
                while(True):
                    line = self.ReadCleanLine()
                    if line is None:
                        break

                    l = line.split()
                    key = l[0]
                    if key == 'DOPTPRM' : continue
                    if key == 'DTPL' : continue
                    if key == 'DCONSTR' : continue
                    if key == 'DCONADD' : continue
                    if key == 'DRESP1' : continue
                    if line[0:5] == 'PARAM' : continue
                    if key == 'CORD2R' : continue
                    if key == '+' : continue
                    if key == 'ENDATA' : break
                    if key == 'PSOLID' : continue
                    if key == 'MAT1' : continue

                    if l[0] == 'SPC' :
                        tag = res.GetNodalTag("SPC"+str(int(line[8:8*2]) ))
                        tag.AddToTag(filetointernalid[int(line[8*2:8*3]) ])
                        continue
                    if key == 'FORCE' :
                        name = "FORCE"+str(int(line[8:8*2]) )
                        tag = res.GetNodalTag(name)
                        node = filetointernalid[int(line[8*2:8*3]) ]
                        tag.AddToTag(node)

                        fx = ParseFloat(line[8*4:8*5])
                        fy = ParseFloat(line[8*5:8*6])
                        fz = ParseFloat(line[8*6:8*7])
                        val = [fx,fy,fz]
                        resdic[name] = {'N':node, "factor":ParseFloat(line[8*3:8*4]), "val":val }
                        continue
                    if key == 'GRID' :
                      line += "             "
                      oid = int(line[8:8*2])
                      x = ParseFloat(line[8*3:8*4])
                      y = ParseFloat(line[8*4:8*5])
                      z = ParseFloat(line[8*5:8*6])
                      ids.append(oid)
                      xs.append(x)
                      ys.append(y)
                      zs.append(z)
                      filetointernalid[oid] = cpt
                      cpt +=1

                      continue
                    if key ==  'CTETRA':
                      oid = int(line[8:8*2])
                      idd2 = str(int(line[8*2:8*3]))
                      P1 = int(line[8*3:8*4])
                      P2 = int(line[8*4:8*5])
                      P3 = int(line[8*5:8*6])
                      P4 = int(line[8*6:8*7])

                      conn = [filetointernalid[xx] for xx in  [P1, P2,P3,P4 ] ]
                      elements = res.GetElementsOfType(ElementNames.Tetrahedron_4)

                      localid = elements.AddNewElement(conn,oid)
                      elements.tags.CreateTag("ET"+idd2,False).AddToTag(localid-1)
                      #res.AddElementToTagUsingOriginalId(oid,)

                      continue

                    print("string '" + str(line) + "' not treated")


            if key == 'ENDATA' : break
            print("string '" + str(line) + "' not treated")
            if line is None: break
            raise

        res.nodes = np.array([xs,ys,zs],dtype=np.float).T
        res.originalIDNodes = np.array(ids,dtype=np.int)
        self.EndReading()
        res.PrepareForOutput()
        return res,resdic

def ParseFloat(line):
    if line[6] == "-"  and line[5] != " ":
        #print(line)
        line = line[0:6] + "e" + line[6:8]
    elif line[5] == "-"  and line[4] != " ":
        #print(line)
        line = line[0:5] + "e" + line[5:8]
    elif line[4] == "-"  and line[3] != " ":
        #print(line)
        line = line[0:4] + "e" + line[4:8]
    elif line[3] == "-"  and line[2] != " ":
        #print(line)
        line = line[0:3] + "e" + line[3:8]

    return float(line)



def CheckIntegrity(GUI = False):


    from BasicTools.Helpers.Tests import TestTempDir

    testData=u"""
BEGIN BULK

GRID      162969        53.00389234.456633.29904
GRID      156839        50.53256233.782536.31633
GRID      156554        50.67131231.871634.13537
GRID      146810        54.21149233.586335.82001
CTETRA    581279       2  162969  156839  156554  146810
SPC            1  162969  1234560.0
FORCE          2  156839       01.0     2949.8  -5792.120.0
ENDDATA
"""


    FR = FemReader()
    FR.SetStringToRead(testData)
    res = FR.Read()[0]

    print(res.nodes)
    print(res)
    from BasicTools.IO.XdmfWriter import WriteMeshToXdmf
    WriteMeshToXdmf(TestTempDir().GetTempPath()+"FemReaderTest.xdmf",res)
    print(TestTempDir().GetTempPath())

    if GUI:
        import os
        os.system('vglrun paraview ' + TestTempDir().GetTempPath()+"FemReaderTest.xdmf")
    return 'ok'

if __name__ == '__main__':
    print(CheckIntegrity(GUI = True))# pragma: no cover
