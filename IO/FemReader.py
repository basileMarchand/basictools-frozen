# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-

import numpy as np
import BasicTools.FE.UnstructuredMesh as UM
from BasicTools.IO.ReaderBase import ReaderBase
import BasicTools.FE.ElementNames as ElementNames


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

            if line[0:8] == 'DESGLB  ' : continue
            if line[0:8] == 'SUBCASE ' : continue
            if line[0:8] == 'ANALYSIS' : continue
            if line[0:8] == '  SPC = ' : continue
            if line[0:8] == '  LOAD =' : continue
            if line[0:8] == 'DESOBJ(M' : continue
            if line[0:8] == '  LABEL ' : continue
            if line == "BEGIN BULK" :
                while(True):
                    line = self.ReadCleanLine()
                    if line[0:8] == 'DOPTPRM ' : continue
                    if line[0:8] == 'DTPL    ' : continue
                    if line[0:8] == 'DCONSTR ' : continue
                    if line[0:8] == 'DCONADD ' : continue
                    if line[0:8] == 'DRESP1  ' : continue
                    if line[0:8] == 'PARAM,CH' : continue
                    if line[0:8] == 'CORD2R  ' : continue
                    if line[0:8] == '+       ' : continue
                    if line[0:6] == 'ENDATA' : break
                    if line[0:8] == 'PSOLID  ' : continue
                    if line[0:8] == 'MAT1    ' : continue
                    if line[0:8] == 'SPC     ' :
                        tag = res.GetNodalTag("SPC"+str(int(line[8:8*2]) ))
                        tag.AddToTag(filetointernalid[int(line[8*2:8*3]) ])
                        continue
                    if line[0:8] == 'FORCE   ' :
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
                    if line[0:8] == 'GRID    ' :
                      #print(line)
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
                    if line[0:8] ==  'CTETRA  ':
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
                    break
                break
            print("string '" + str(line) + "' not treated")
            break

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

def ReadFem(fileName=None,string=None):
    obj = FemReader()
    if fileName is not None:
        obj.SetFileName(fileName)

    if string is not None:
        obj.SetStringToRead(string)

    return obj.Read()[0]

def CheckIntegrity(GUI = False):


    from BasicTools.Helpers.Tests import TestTempDir
    newFileName = "/home/fbordeu-weld/PythonTools/cas_application_ASL/cas_application_ASL.fem"

    res = ReadFem(fileName = newFileName)
    print(res.nodes)
    print(res)
    from BasicTools.IO.XdmfWriter import WriteMeshToXdmf
    WriteMeshToXdmf(TestTempDir().GetTempPath()+"FemReaderTest.xdmf",res)
    print(TestTempDir().GetTempPath())

    if GUI:
        import os
        os.system('vglrun /home/software/paraview/binv5.2.0/bin/paraview  ' + TestTempDir().GetTempPath()+"FemReaderTest.xdmf")

    return 'ok'

if __name__ == '__main__':
    print(CheckIntegrity(GUI = True))# pragma: no cover
