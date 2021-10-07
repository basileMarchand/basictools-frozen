# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
                       

""" Inp file reader (Abaqus simulation file)
"""
import numpy as np
import csv

from BasicTools.Helpers.Timer import Timer
import BasicTools.Containers.ElementNames as EN
from BasicTools.Containers.Filters import ElementFilter
from BasicTools.IO.ReaderBase import ReaderBase
from BasicTools.IO.AbaqusTools import InpNameToBasicTools, permutation

KeywordToIgnore = ["INITIAL CONDITIONS",
                   "AMPLITUDE",
                   "EXPANSION",
                   "DISTRIBUTION TABLE",
                   "COUPLING",
                   "SOLID SECTION",
                   "CONNECTOR SECTION",
                   "SURFACE",
                   "SURFACE INTERACTION",
                   "FRICTION",
                   "CONTACT PAIR",
                   "CLEARANCE",
                   "PARAMETER",
                   "PART",
                   "END PART",
                   "ASSEMBLY",
                   "INSTANCE",
                   "END INSTANCE",
                   "END ASSEMBLY"
                   ]

def JoinInp(filename):
    fII = open("join_"+filename,"w")

    def copycontent(filename):
        for l in open(filename,"r"):
            if len(l) >=8 and l[0:8]  == "*INCLUDE":
                copycontent( l.split("=")[1].strip()) 
                pass
            else:
                fII.write(l)

    copycontent(filename)
    
def SplitInp(filename):
    import os
    f = InpReader()
    f.SetFileName(fileName= filename)
    fII = open("split_"+filename,"w")
    f.StartReading()
    paires = {"PART":"*END PART",
              "INSTANCE":"*END INSTANCE",
              "STEP":"*END STEP",
              "Step":"*End Step",
              "ASSEMBLY":"*END ASSEMBLY",}


    def writeonfile(inputfile,outputfile,waitfor="*",pwd="./"):
        print("working on " , " ", pwd)
        if pwd != "./" and waitfor != "*":
            os.mkdir(pwd)
        kcpt = 0
        l = inputfile.ReadCleanLine()
        while(True):
            if l is None:
                return None

            if l[0:len(waitfor)] == waitfor:
                if waitfor == "*":
                    return l
                else:
                    outputfile.write(l+"\n")
                    #print("geting out",waitfor)
                    return  inputfile.ReadCleanLine()

            if l[0] == "*":
                ldata = LineToDic(l)
                keyword = ldata["KEYWORD"]+"_"+str(kcpt)
                kcpt += 1

                outputfile.write("*INCLUDE INPUT="+pwd+keyword+".partial\n")
                kf = open(pwd+keyword+".partial","w")

                kf.write(l+"\n")
                if ldata["KEYWORD"] in paires:
                    waitforin=paires[ldata["KEYWORD"]]
                else:
                    waitforin="*"
                print("Creating file ", pwd+keyword+"/main.partial","w")

                l = writeonfile(f,kf,waitfor=waitforin,pwd=pwd+keyword+"/")
            else:
                outputfile.write(l+"\n")
                l = inputfile.ReadCleanLine()

    writeonfile(f,fII,waitfor="$$",pwd="./")

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
                  res[s[0].lstrip().rstrip().upper()] = s[1].lstrip().rstrip().lstrip('"').rstrip('"')
                else:
                  res[f.lstrip().rstrip().upper()] = True
    return res

def LineToList(text):
   return  list(csv.reader([text], delimiter=',', quotechar='"'))[0]

def LineToListNoQuote(text):
    return [s.strip() for s in text.split(",")]

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

    def find(self,l,expr):
        res = l.upper().find(expr.upper())
        return res

    def ReadCleanLine(self):
        # in the inp abaqus the lines ending in ',' are splited in multiple lines
        res =  super(InpReader,self).ReadCleanLine()
        if res is None:
            return res
        if res[-1] == ",":
            nline = self.PeekLine()
            if nline[0] != "*":
                return res + " "+  super(InpReader,self).ReadCleanLine()
        return res 

    def Read(self,fileName=None,string=None, out=None):
        import BasicTools.FE.ProblemData as ProblemData
        from BasicTools.Linalg.Transform  import Transform
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

            if  "KEYWORD" in ldata and ldata["KEYWORD"] in KeywordToIgnore:
                print("Ignoring keyword: '" + str(ldata["KEYWORD"]) +"' jumping to the next star" )
                l = DischardTillNextStar(self.ReadCleanLine )
                continue

            if self.find(l,"*ELSET")>-1:
                elsetName = ldata['ELSET']
                l  = self.ReadCleanLine()
                ds = []
                tagsNames = []
                while(True):
                    if l is None or l.find("*") > -1 or not l:
                        tags = {}
                        for d in ds:
                            elems, id = filetointernalidElement[d]
                            if elems.elementType not in tags:
                                tag = elems.tags.CreateTag(elsetName)
                                tags[elems.elementType] = tag
                            else:
                                tag = tags[elems.elementType]
                            tag.AddToTag(id)
                        if len(tagsNames) >0 :
                            for elname,eldata,elid in ElementFilter(res,tags=tagsNames):
                                tag = eldata.tags.CreateTag(elsetName,False).AddToTag(elid)
                        break
                    s = LineToListNoQuote(l)
                    for ss in s:
                        try:
                            ds.append(int(ss))
                        except:
                            tagsNames.append(ss)
                    l = self.ReadCleanLine()
                    continue
                continue

            if self.find(l,"**LENGTH UNITS")>-1:
                if l.find("mm")>-1:
                    coef = 0.001
                l = self.ReadCleanLine()
                continue

            if self.find(l,"*NODE")>-1:
                    nodes= []
                    originalIds = []

                    cpt = 0
                    res.originalIDNodes = np.empty((0,1), int)
                    res.nodes = np.empty((0,3), float)
                    l = self.ReadCleanLine()
                    dim = 3
                    s = None

                    while(True):
                        if len(l) == 0 or l[0]=="*" or not l:
                            break
                        s = LineToListNoQuote(l)
                        oid = int(s[0])
                        filetointernalid[oid] = cpt
                        nodes.append(list(map(float,s[1:])) )
                        originalIds.append(oid)
                        cpt += 1
                        l = self.ReadCleanLine()

                    if s is not None:
                        # we use the last point to detect the dimensionality of the mesh 
                        dim = len(s)-1
                    res.nodes = np.array(nodes,dtype=float)
                    res.nodes.shape = (cpt,dim)
                    res.originalIDNodes = np.array(originalIds,dtype=int)
                    continue

            if self.find(l,"*ELEMENT")>-1:
                data = LineToDic(l)
                etype = data["TYPE"]
                nametype = InpNameToBasicTools[etype]
                per = permutation.get(etype,None)

                l = self.ReadCleanLine()
                elements = res.GetElementsOfType(nametype)
                initialcid = elements.GetNumberOfElements()

                while(True):
                    if l is None or not l or l[0]== "*" :
                        break
                    s = LineToListNoQuote(l)
                    oid = int(s[0])
                    conn = [filetointernalid[x] for x in  map(int,s[1:]) ]

                    cid = elements.AddNewElement(conn,oid)-1
                    filetointernalidElement[oid] = (elements,cid)

                    l = self.ReadCleanLine()
                    
                if per is not None:
                    elements.connectivity = elements.connectivity[:,per]

                finalcid = elements.GetNumberOfElements()

                if "ELSET" in data:
                    elset = data["ELSET"]
                    elements.GetTag(elset).AddToTag(range(initialcid,finalcid))
                elif etype == "CONN3D2":
                    elements.GetTag("CONN3D2").AddToTag(range(initialcid,finalcid))

                res.ComputeGlobalOffset()
                continue
         
            if self.find(l,"*NSET")>-1:
                data = LineToDic(l)
                nsetName = data['NSET']
                l  = self.ReadCleanLine()
                nset = []

                if "GENERATE" in ldata and ldata["GENERATE"]:
                    d = LineToList(l)
                    d = (list(map(int,d) ))
                    nset = range(d[0],d[1]+1,d[2])
                    tag = res.nodesTags.CreateTag(nsetName,False)    
                    tagsids = [filetointernalid[x] for x in  nset if x in filetointernalid ]
                    if len(tagsids) != len(nset):
                        print("Warning NSET GENERATE : not all the nset generated are present in the mesh ")
                    tag.AddToTag(tagsids) 
                    l = self.ReadCleanLine()
                    continue
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

                tag = res.nodesTags.CreateTag(nsetName,False)
                tag.AddToTag([filetointernalid[x] for x in  nset ])
                continue

            if self.find(l,"*DISTRIBUTION")>-1:
                data = LineToDic(l)
                location = data['LOCATION'] # "ELEMENT"
                name = data['NAME'] # "ELEMENT"
                table = data['TABLE'] # "ELEMENT"

                if table == "DistributionTable_Orientation":
                    NumberOfData = 6
                elif table == "DistributionTable_Density":
                    NumberOfData = 1
                elif table == "DistributionTable_Elastic":
                    NumberOfData = 9
                elif table == "DISTRIB_TABLE":
                    NumberOfData = 9
                else:
                    raise Exception("Error! ")
                data = np.zeros((res.GetNumberOfElements(),NumberOfData))
                if location == "ELEMENT":
                    l = self.ReadCleanLine()
                    s = l.split(',')
                    if len(s[0]) == 0:
                        #default values
                        s = list(map(float,s[1:]))
                        if len(s) < NumberOfData:
                            # multiline data
                            s.extend(map(float,self.ReadCleanLine().split(',')) )
                        s = np.array(s)
                        for i,v in enumerate(s):
                            data[:,i] = v

                    l = self.ReadCleanLine()
                    while(True):
                        if l is None or l.find("*") > -1 or not l:
                            break

                        s = l.split(",")
                        elset = s[0].split(".")[-1]
                        s = list(map(float,s[1:]))
                        if len(s) < NumberOfData:
                            # multiline data
                            s.extend(map(float,self.ReadCleanLine().replace(',', ' ').split()) )

                        s= np.array(s)
                        for elname,eldata,elids in ElementFilter(res,tag=elset):
                            data[elids+eldata.globaloffset,:] = s

                        l = self.ReadCleanLine()
                    
                    if table == "DistributionTable_Orientation":
                        res.elemFields["V1"] = data[:,0:3]
                        res.elemFields["V2"] = data[:,3:]
                    elif table == "DistributionTable_Density":
                        res.elemFields["density"] = data
                    elif table == "DistributionTable_Elastic":
                        res.elemFields["E"] = data
                        res.elemFields["E1"]   = data[:,0]
                        res.elemFields["E2"]   = data[:,1]
                        res.elemFields["E3"]   = data[:,2]
                        res.elemFields["Nu12"] = data[:,3]
                        res.elemFields["Nu13"] = data[:,4]
                        res.elemFields["Nu23"] = data[:,5]
                        res.elemFields["G12"]  = data[:,6]
                        res.elemFields["G13"]  = data[:,7]
                        res.elemFields["G23"]  = data[:,8]
                    
                continue
              
            if self.find(l,"*HEADING")>-1:
                l = self.ReadCleanLine()
                HEADING = ""
                while(True):
                    if l is None or l.find("*") > -1 or not l:
                        break
                    HEADING += str(l) + "\n"
                    l = self.ReadCleanLine()

                continue

            if self.find(l,"*MATERIAL") > -1:
                data = LineToDic(l)
                name = data["NAME"]
                mat = ProblemData.Material()
                while(True):
                    l = self.ReadCleanLine()
                    if l is None:
                        break
                    data = LineToDic(l)
                    t = data["KEYWORD"]
                    if t not in ["ELASTIC","DENSITY","EXPANSION","ELASTIC","DISTRIBUTION"] :
                        break

                    if t == "DISTRIBUTION":
                        print("Ignoring DISTRIBUTION")
                        l = DischardTillNextStar(self.ReadCleanLine )
                        break

                    if t == "EXPANSION":
                        if "ZERO" in data: 
                            mat.AddProperty(t,"ZERO",data["ZERO"])

                    if "TYPE" in data:
                        st = data["TYPE"]
                    else:
                        st = None
                    l = self.ReadCleanLine()
                    s = list(l.strip('\n').strip().split() )
                    mat.AddProperty(t,st,s)

                meta.materials[name] = mat
                continue

            if self.find(l,"*ORIENTATION")>-1:
                data = LineToDic(l)
                orient = Transform()

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

            if self.find(l,"*SURFACE")>-1:
                data = LineToDic(l)
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
                        elements,cid = filetointernalidElement[originalElemNumber][0]
                        connectivity = elements.connectivity[cid,:]

                        typeAndConnectivity = EN.faces[elements.elementType][faceNumber]
                        faceconn = connectivity[typeAndConnectivity[1]]

                        elements = res.GetElementsOfType(typeAndConnectivity[0])
                        cid = elements.AddNewElement(faceconn,-1)

                        elements.GetTag(name).AddToTag(cid-1)
                        l = self.ReadCleanLine()

                        continue
                    res.ComputeGlobalOffset()
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


            print(l)
            print(ldata)
            print(("line starting with <<"+l[:20]+">> not considered in the reader"))
            l = self.ReadCleanLine()
            raise Exception("Error reading file")

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
