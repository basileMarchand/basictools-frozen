# -*- coding: utf-8 -*-
""" Fem file reader

"""

import numpy as np

__author__ = "Felipe Bordeu"
import BasicTools.Containers.UnstructuredMesh as UM
import BasicTools.Containers.ElementNames as ElementNames

from BasicTools.IO.ReaderBase import ReaderBase


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
        self.SetFileName(fileName)
        self.ignoreNotTreated = True

    def GetField(self,line,number):
        return line[8*number:8*(number+1)]

    def Read(self):
        self.StartReading()

        res = UM.UnstructuredMesh()
        resdic = {}
        resdic["CordinateSystems"] = {}
        resdic["Objective"] = None
        resdic["Optim Functions"] = {}
        resdic["Cases"] = {}
        ids = []
        xs =[]
        ys = []
        zs = []
        filetointernalid = {}
        cpt =0

        elements = res.GetElementsOfType(ElementNames.Tetrahedron_4)
        elements.Reserve(1000000)
        need_to_read = True
        while(True):
            if need_to_read :
                line = self.ReadCleanLine()
                need_to_read = True
            if line is None:
                break

            l = line.split()
            key = l[0]
            #print(key)
            if key == 'DESGLB' and self.ignoreNotTreated :  continue
            if key == 'SUBCASE' :
            #https://knowledge.autodesk.com/support/nastran/learn-explore/caas/CloudHelp/cloudhelp/2019/ENU/NSTRN-Reference/files/GUID-38938C32-4019-4DDA-AE94-2987C8355669-htm.html
                subcase = int(self.GetField(line,1))
                data = {}

                while(True):
                  line = self.ReadCleanLine()
                  key = line.split()[0]
                  #print("tt" +line)
                  if line.find("LABEL") > -1:
                  #https://knowledge.autodesk.com/support/nastran/learn-explore/caas/CloudHelp/cloudhelp/2019/ENU/NSTRN-Reference/files/GUID-99CAEA1D-352E-4CDA-B74B-C58E4D12D020-htm.html?st=label
                     data['label'] = line.split("LABEL")[1].strip()
                     continue
                  if key == 'ANALYSIS':
                     #https://knowledge.autodesk.com/support/nastran/learn-explore/caas/CloudHelp/cloudhelp/2019/ENU/NSTRN-Reference/files/GUID-E7589B2F-C046-4B34-A668-B4E796172645-htm.html
                     #print(line)
                     data['analisys'] = line.split()[1].strip()
                     continue
                  if key == 'SPC' :
                     #https://knowledge.autodesk.com/support/nastran/learn-explore/caas/CloudHelp/cloudhelp/2019/ENU/NSTRN-Reference/files/GUID-A456B121-CEB2-48B6-9EBA-7F19EF250AA5-htm.html
                     data['SPC'] = int(line.split("=")[1].strip())
                     continue
                  if key == 'LOAD':
                     # External Static Load Set Selection
                     #https://knowledge.autodesk.com/support/nastran/learn-explore/caas/CloudHelp/cloudhelp/2019/ENU/NSTRN-Reference/files/GUID-02DF840A-A5EB-4F04-AB79-DC82A5691105-htm.html
                     data['LOAD'] = int(line.split("=")[1].strip())
                     continue

                  if key == 'MPCFORCE'  :
                      print("MPCFORCE intensionally ignored")
                      continue
                     # OUTPOUT
                     #  Requests multipoint constraint force vector output
                     #https://knowledge.autodesk.com/support/nastran/learn-explore/caas/CloudHelp/cloudhelp/2019/ENU/NSTRN-Reference/files/GUID-97B55735-D7B3-4396-A161-BC6B8B694B98-htm.html?st=MPCFORCES

                  if key == 'STRAIN'  :
                      print("STRAIN intensionally ignored")
                      continue
                     # OUTPOUT
                     #  Element Strain Output Request
                     #https://knowledge.autodesk.com/support/nastran/learn-explore/caas/CloudHelp/cloudhelp/2019/ENU/NSTRN-Reference/files/GUID-5DE89858-947A-4B35-B244-03BEE2CA779C-htm.html
                  if key == 'DESSUB'  :
                     # Subcase Level Design Constraints Selection
                     #https://knowledge.autodesk.com/support/nastran/learn-explore/caas/CloudHelp/cloudhelp/2019/ENU/NSTRN-Reference/files/GUID-093E2B66-D5C4-4C2E-8494-7E6AD36F8B01-htm.html
                     data['DESSUB'] = int(line.split("=")[1].strip())
                     continue
                  break
                resdic['Cases'][subcase] = data
                need_to_read = False
                continue

            if line[0:6] == 'DESOBJ' :
                vec = line.split("=")
                resdic["Objective"] = {"id":int(vec[1]),"type":vec[0].split("(")[1].split(")")[0]}
                continue
            if key == 'LABEL' and self.ignoreNotTreated : continue

            if key == "SET":
                need_to_read = False
                l = line.split()
                key = l[0]
                print("--> " + key)
                setNumber = l[1]

                if l[2] != "=":
                    raise

                line = "".join(l[3:])

                while(True):
                    if need_to_read :
                        line = self.ReadCleanLine()
                        need_to_read = True
                    if line is None:
                        break

                    units = line.split(",")
                    for u in units:
                        if len(u) == 0:
                            continue
                        if "THRU" in u:
                            partition = u.split()
                            begin = partition[0]
                            end = partition[-1]
                        else:
                            idd = int(u)
                    need_to_read = True
                    if line[-1] != ",":
                        break
                continue


            if line == "BEGIN BULK" :
                need_to_read = True
                while(True):
                    if need_to_read :
                        line = self.ReadCleanLine()
                        need_to_read = True
                    if line is None:
                        break

                    l = line.split()
                    key = l[0]
                    #print("---> " + key)
                    if key == 'DOPTPRM'  :
                        resdic["DOPTPRM"] = {self.GetField(line,1):ParseFloat(self.GetField(line,2))}
                        continue
                    if key == 'DTPL' :
                        print("DTPL intensionally ignored")
                        continue
                    if key == 'DCONSTR' and self.ignoreNotTreated :
                        need_to_read    = True
                        continue
                    if key == 'DCONADD' and self.ignoreNotTreated : continue
                    if key == 'DRESP1' :
                        #https://knowledge.autodesk.com/support/nastran/learn-explore/caas/CloudHelp/cloudhelp/2019/ENU/NSTRN-Reference/files/GUID-4C44D4BF-E70D-4548-8BED-D4C0497E5479-htm.html?st=DRESP1
                        idd = int(self.GetField(line,1))
                        userlabel = self.GetField(line,2).strip()
                        RTYPE = self.GetField(line,3).strip()
                        REGION = self.GetField(line,5).strip()
                        ATTA = self.GetField(line,6).strip()
                        ATTB = self.GetField(line,7).strip()
                        ATTC = self.GetField(line,8).strip()
                        bulkdata = []
                        line = self.ReadCleanLine()
                        while line[0] == "+":
                            try:
                                #print("--" + line)
                                for x in [1,2,3,4,5,6,7,8]:
                                   bulkdata.append(int(self.GetField(line,x)))
                                line = self.ReadCleanLine()
                            except Exception as inst:
                                #print(inst)
                                line = self.ReadCleanLine()
                                break
                        need_to_read    = False
                        #print(bulkdata)
                        continue

                    if line[0:5] == 'PARAM' and self.ignoreNotTreated : continue
                    if key == 'CORD2R' :
                        # from
                        # https://knowledge.autodesk.com/support/nastran/learn-explore/caas/CloudHelp/cloudhelp/2018/ENU/NSTRN-Reference/files/GUID-F048B85B-79B9-49AD-94B4-87CC69122C3B-htm.html
                        idd = int(self.GetField(line,1))
                        floats = [ParseFloat(self.GetField(line,x)) for x in [3,4,5,6,7,8]]
                        line = self.ReadCleanLine()
                        floats.extend([ParseFloat(self.GetField(line,x)) for x in [1,2,3]] )
                        from BasicTools.FE.ProblemData import Orientation
                        orient = Orientation()
                        orient.system = "CORD2R"

                        A = np.array(floats[0:3])
                        B = np.array(floats[3:6])
                        C = np.array(floats[6:9])
                        z = (B-A)/np.linalg.norm(B-A)
                        x = (C-A)- z*np.dot(z,C-A)
                        y = np.cross(z,x)

                        orient.SetOffset(A)
                        orient.SetFirst(x)
                        orient.SetSecond(y)
                        resdic["CordinateSystems"][idd] = orient

                        continue
                    if key == 'CORD1R' and self.ignoreNotTreated : continue
                    if key == 'RBE2' and self.ignoreNotTreated : continue
                    if key == 'RBE3' and self.ignoreNotTreated : continue

                    if key == '+'  :
                        #print(line )
                        print("+ intensionally ignored")
                        continue
                    if key == 'ENDATA' : break
                    if key == 'PSOLID' and self.ignoreNotTreated : continue
                    if key == 'MAT1' and self.ignoreNotTreated : continue
                    if key == 'LOADADD' and self.ignoreNotTreated : continue

                    if l[0] == 'SPC' :
                        tag = res.GetNodalTag("SPC"+str(int(line[8:8*2]) ))
                        tag.AddToTag(filetointernalid[int(line[8*2:8*3]) ])
                        continue
                    if key == 'PLOAD4' :
                        name = "PLOAD4"+str(int(line[8:8*2]) )
                        tag = res.GetNodalTag(name)
                        node = filetointernalid[int(line[8*2:8*3]) ]
                        tag.AddToTag(node)
                        resdic[name] = {'N':node, "val":ParseFloat(line[8*3:8*4]) }
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

                      oid = int(line[8:8*2])

                      x = ParseFloat(line[8*3:8*4])
                      y = ParseFloat(line[8*4:8*5])
                      z = ParseFloat(line[8*5:8*6])
                      ids.append(oid)
                      csystem= self.GetField(line,2)
                      if csystem != " "*8:
                          x,y,z = resdic["CordinateSystems"][int(csystem)].GetPointInOriginalSystem([x,y,z])
                      xs.append(x)
                      ys.append(y)
                      zs.append(z)
                      filetointernalid[oid] = cpt
                      cpt +=1

                      continue
                    if key ==  'CTRIA3':
                      oid = int(line[8:8*2])
                      idd2 = str(int(line[8*2:8*3]))
                      P1 = int(line[8*3:8*4])
                      P2 = int(line[8*4:8*5])
                      P3 = int(line[8*5:8*6])

                      conn = [filetointernalid[xx] for xx in  [P1, P2,P3 ] ]
                      elements = res.GetElementsOfType(ElementNames.Triangle_3)

                      localid = elements.AddNewElement(conn,oid)
                      elements.tags.CreateTag("ET"+idd2,False).AddToTag(localid-1)
                      #res.AddElementToTagUsingOriginalId(oid,)
                      continue
                    if key ==  'CTETRA':
                      oid = int(line[8:8*2])
                      idd2 = str(int(line[8*2:8*3]))
                      points = []
                      for i in range(3,9):
                          try:
                             points.append(int(line[8*i:8*(i+1)]))
                          except:
                             pass

                      if len(points) > 4:
                          line = self.ReadCleanLine()
                          for i in range(1,5):
                             points.append(int(line[8*i:8*(i+1)]))

                      conn = [filetointernalid[xx] for xx in  points ]

                      if len(conn) == 4:
                          elements = res.GetElementsOfType(ElementNames.Tetrahedron_4)
                      elif len(conn) == 10:
                          elements = res.GetElementsOfType(ElementNames.Tetrahedron_10)

                      localid = elements.AddNewElement(conn,oid)
                      elements.tags.CreateTag("ET"+idd2,False).AddToTag(localid-1)
                      #res.AddElementToTagUsingOriginalId(oid,)
                      continue
                    if key == 'ENDDATA' : break
                    if line is None: break

                    print("string '" + str(line) + "' not treated")
                    raise(ValueError("string  '" + str(line) + "' not treated"))

            if key == 'ENDDATA' : break
            if line is None: break
            raise(ValueError("string '" + str(line) + "' not treated"))

        res.nodes = np.array([xs,ys,zs],dtype=np.float).T
        res.originalIDNodes = np.array(ids,dtype=np.int)
        self.EndReading()
        res.PrepareForOutput()
        return res,resdic

from BasicTools.IO.IOFactory import RegisterReaderClass
RegisterReaderClass(".fem",FemReader)


def ParseFloat(line):
    line += (8-len(line))*' '
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
