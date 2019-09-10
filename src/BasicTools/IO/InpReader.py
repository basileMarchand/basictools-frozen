# -*- coding: utf-8 -*-

""" Inp file reader (Abaqus simulation file)
"""
import numpy as np

import BasicTools.Containers.ElementNames as EN
from BasicTools.IO.ReaderBase import ReaderBase


KeywordToIgnore = ["INITIAL CONDITIONS",
                   "AMPLITUDE",
                   "EXPANSION",
                   "DISTRIBUTION TABLE",
                   "DISTRIBUTION",
                   "COUPLING",
                   "SOLID SECTION",
                   "CONNECTOR SECTION"
                   ]

InpNumber = {}
InpNumber['C3D4'] = EN.Tetrahedron_4
InpNumber['S3'] = EN.Triangle_3
InpNumber['CONN3D2'] = EN.Bar_2

def LineToDic(text):
    import csv
    res = {}
    for l in csv.reader([text], delimiter=',', quotechar='"'):
        for f in l:
            if len(f) == 0:
                continue
            if f[0] == "*":
                res["KEYWORD"] = f[1:]
            else:
                if f.find("=") >-1:
                  s = f.split("=")
                  res[s[0].lstrip().rstrip()] = s[1].lstrip().rstrip().lstrip('"').rstrip('"')
                else:
                  res[f.lstrip().rstrip()] = True
    return res

def LineToList(text):
   import csv
   return  list(csv.reader([text], delimiter=',', quotechar='"'))[0]

def ReadLine(string):
  raise
  while(True):
     l = string.readline().strip('\n').lstrip().rstrip()
     if len(l) >= 2:
        ##comment
        if l[0:2] == "**":
           continue
     return l

def DischardTillNextStar(func):
    while(True):
        courrentText = func()
        if courrentText is None:
            break
        if len(courrentText) > 1 and courrentText[0] == "*":
           break

    return  courrentText

def ReadInp(fileName=None,string=None,out=None,**kwargs):
    reader = InpReader()
    reader.SetFileName(fileName)
    reader.SetStringToRead(string)
    reader.Read(fileName=fileName, string=string,out=out,**kwargs)
    return reader.output

class InpReader(ReaderBase):
    def __init__(self):
        super(InpReader,self).__init__()
        self.commentChar= "**"
        self.readFormat = 'r'

    def Read(self,fileName=None,string=None, out=None):
        import BasicTools.FE.ProblemData as ProblemData
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

        coef = 1.
        meta = ProblemData.ProblemData()

        filetointernalid = {}
        filetointernalidElement = {}
        l = self.ReadCleanLine()

        while(True):
          print(l)
          if not l: break
          ldata = LineToDic(l)


          if l.find("**LENGTH UNITS")>-1:
              if l.find("mm")>-1:
                  coef = 0.001
              l = self.ReadCleanLine()
              continue

          if l.find("*NODE")>-1:
                nodes= []
                originalIds = []

                cpt = 0
                res.originalIDNodes = np.empty((0,1), int)
                res.nodes = np.empty((0,3), float)
                l = self.ReadCleanLine()
                while(True):
                    if len(l) == 0:
                        continue
                    if l.find("*") > -1 or not l:
                        break
                    s = l.replace(',', '').split()
                    oid = int(s[0])
                    if oid not in filetointernalid:
                      filetointernalid[oid] = cpt
                      #res.originalIDNodes = np.vstack((res.originalIDNodes,int(s[0])))
                      #res.nodes = np.vstack((res.nodes,list(map(float,s[1:]))))
                      nodes.append(list(map(float,s[1:])) )
                      originalIds.append(int(s[0]) )
                      cpt += 1

                    l = self.ReadCleanLine()

                res.nodes = np.array(nodes,dtype=np.float)
                res.nodes.shape = (cpt,3)
                res.originalIDNodes = np.array(originalIds,dtype=np.int)

                continue

          if l.find("*ELEMENT")>-1:
            data = LineToDic(l)
            #s = l.replace(',', '').split()
            etype = data["TYPE"]
            nametype = InpNumber[etype]
            hasElset = False
            if "ELSET" in data:
                elset = data["ELSET"]
                hasElset = True
            l = self.ReadCleanLine()
            while(True):
              if l is None:
                  break
              if l.find("*") > -1 or not l:
                   break
              s = l.replace(',', ' ').split()
              conn = [filetointernalid[x] for x in  map(int,s[1:]) ]
              elements = res.GetElementsOfType(nametype)
              oid = int(s[0])
              cid = elements.AddNewElement(conn,oid)
              filetointernalidElement[oid] = (elements,cid-1)
              if hasElset:
                  elements.GetTag(elset).AddToTag(cid-1)

              if etype == "CONN3D2":
                  elements.GetTag("CONN3D2").AddToTag(cid-1)

              #cpt += 1
              l = self.ReadCleanLine()
            continue




          if l.find("*NSET")>-1:
            data = LineToDic(l)
            nsetName = data['NSET']
            l  = self.ReadCleanLine()
            nset = []

            if "GENERATE" in ldata and ldata["GENERATE"]:
                d = LineToList(l)
                d = (list(map(int,d) ))
                nset = range(d[0],d[1]+1,d[2])
                l = self.ReadCleanLine()
            else:
                while(True):
                    if l is None:
                        break
                    if l.find("*") > -1 or not l:
                        break
                    s = l.replace(',', '').split()
                    nset.extend(map(int,s))

                    l = self.ReadCleanLine()
                    continue

            tag = res.nodesTags.CreateTag(nsetName)
            tag.SetIds([filetointernalid[x] for x in  nset ])
            continue

          if l.find("*ELSET")>-1:
            data = LineToDic(l)
            elsetName = data['ELSET']
            l  = self.ReadCleanLine()
            while(True):
                if l is None:
                    break

                if l.find("*") > -1 or not l:
                    break
                s = l.replace(',', '').split()
                for d in map(int,s):
                    #res.AddElementToTagUsingOriginalId(d,elsetName)
                    ElementsAndNumber = filetointernalidElement[d]
                    ElementsAndNumber[0].AddElementToTag(ElementsAndNumber[1],elsetName)
                l = self.ReadCleanLine()
                continue

            continue

          if l.find("*HEADING")>-1:
               HEADING = self.ReadCleanLine()
               meta.HEADING = HEADING
               l = self.ReadCleanLine()
               continue
    #*MATERIAL, NAME="ALU ALsi10mg"
    #*ELASTIC, TYPE=ISOTROPIC
    #75000000000., 0.33
    #*DENSITY
    #2670.

          if l.find("*MATERIAL") > -1:
               data = LineToDic(l)
               name = data["NAME"]
               mat = ProblemData.Material()


               while(True):
                 l = self.ReadCleanLine()
                 data = LineToDic(l)
                 t = data["KEYWORD"]
                 if not (t == "ELASTIC" or t == "DENSITY" or t == "EXPANSION" ) :
                     break

                 if t == "EXPANSION":
                     mat.AddProperty(t,"ZERO",data["ZERO"])

                 if "TYPE" in data:
                    st = data["TYPE"]
                 else:
                     st = None
                 l = self.ReadCleanLine()
                 s = list(map(float, l.strip('\n').lstrip().rstrip().replace(","," ").split() ))
                 mat.AddProperty(t,st,s)

               meta.materials[name] = mat
               continue
    #line to delete
    #           break

          if l.find("*ORIENTATION")>-1:
               data = LineToDic(l)
               orient = ProblemData.Transform()

               name = data["NAME"]
               orient.system = data["SYSTEM"]
               if orient.system != "RECTANGULAR":
                   s = list(map(float, self.ReadCleanLine().replace(","," ").split() ))
                   orient.SetFirst(s[0:3])
                   orient.SetSecond(s[3:6])
                   if len(s) >= 9:
                       orient.SetOffset(s[6:9])

               meta.orientations[name] = orient
               l = DischardTillNextStar(self.ReadCleanLine )
               continue

          if l.find("*SURFACE")>-1:
            data = LineToDic(l)
            #s = l.replace(',', '').split()
            name = data["NAME"]

            if data["TYPE"] == "ELEMENT":
                l = self.ReadCleanLine()

                while(True):
                   if l is None:
                       break
                   if l.find("*") > -1 or not l:
                       break
                   s = l.split(",")
                   originalElemNumber = int(s[0])
                   faceNumber = int(s[1].lstrip().rstrip()[-1])-1
                   elements = filetointernalidElement[originalElemNumber][0]
                   cid = filetointernalidElement[originalElemNumber][1]
                   connectivity = elements.connectivity[cid,:]

                   typeAndConnectivity = EN.faces[elements.elementType][faceNumber]
                   faceconn = connectivity[typeAndConnectivity[1]]

                   elements = res.GetElementsOfType(typeAndConnectivity[0])
                   cid = elements.AddNewElement(faceconn,-1)

                   elements.GetTag(name).AddToTag(cid-1)
                   l = self.ReadCleanLine()

                   continue
                continue
            else:
                raise Exception("NOT IMPLEMENTED sorry")

    #*COUPLING, CONSTRAINT NAME="Coupling-1", REF NODE=239100, SURFACE="Rigid Connection1-1"
    #*KINEMATIC
    #1, 6

          if ldata["KEYWORD"] == "STEP":
             data = LineToDic(l)
             name = data["NAME"]

             l = self.ReadCleanLine()
             cs = ProblemData.StudyCase()

             while(True):
                data = LineToDic(l)


                if "KEYWORD" in data and data["KEYWORD"] == "END STEP":
                   l = self.ReadCleanLine()
                   break

                if "KEYWORD" in data and data["KEYWORD"] == "STATIC":
                   cs.type = "STATIC"

                l = DischardTillNextStar(self.ReadCleanLine)
             continue


          if ldata["KEYWORD"] in KeywordToIgnore:
              l = DischardTillNextStar(self.ReadCleanLine )
              continue
         ##  ----------------------------------------

          print(l)
          print(ldata)
          print(("line starting with <<"+l[:20]+">> not considered in the reader"))
          l = self.ReadCleanLine()
          raise
          continue

        self.EndReading()
        res.originalIDNodes = np.squeeze(res.originalIDNodes)
        res.PrepareForOutput()
        self.output = (res,meta)
        return res


from BasicTools.IO.IOFactory import RegisterReaderClass
RegisterReaderClass(".inp",InpReader)


def CheckIntegrity():
    res1 = LineToDic('*NSET, NSET="Fixed Displacement1", KEY2=5')
    if res1['KEY2'] != '5':
        return "not ok"

    data = u"""** -------------------------------------------------------
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
    print(res)





    return 'ok'

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
