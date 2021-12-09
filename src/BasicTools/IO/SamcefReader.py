# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

import numpy as np

import BasicTools.Containers.ElementNames as EN
from BasicTools.IO.ReaderBase import ReaderBase
from BasicTools.Helpers import ParserHelper as PH
from BasicTools.NumpyDefs import PBasicFloatType, PBasicIndexType

KeywordToIgnore = ["INIT",
                   "MCNL",
                   "ASEF",
                   "HYP",
                   "DES",
                   "MAT",
                   "UNIT",
                   "GEL",
                   "SAM",
                   "OPT",
                   "SAI",
                   "UNITE",
                   "PHP",
                   "RBE",
                   "BPR",
                   "FCT",
                   "SUB",
                   "MCT",
                   ]

def DischardTillNextSection(func):
    while(True):
        courrentText = func()
        if courrentText is None:
            break
        if len(courrentText) > 1 and courrentText[0] == ".":
           break

    return  courrentText

def LineToList(text):
   import csv
   return  list(csv.reader([text], delimiter=' ', quotechar='"'))[0]

def LineToDic(text,res = None):
    if res is None:
        res = {}

    res["&"] = False
    fields = LineToList(text)
    cpt = 0
    if len(fields[0]) > 0  and fields[0][0] == ".":
        res["KEYWORD"] = text.split()[0].split(".")[1]
        cpt +=1
    else:
        res["KEYWORD"] = None

#    ignored = []
    while(cpt < len(fields)):
        k = fields[cpt]
        if k in res.keys():
            cpt +=1
            if type(res[k]) == bool:
                data = True
            else:
                if type(res[k]).__module__ == np.__name__:
                    l = len(res[k])
                    data = PH.Read(fields[cpt:cpt+l],res[k])
                    cpt += l
                else:
                    data = PH.Read(fields[cpt],res[k])
                    cpt +=1
            res[k] = data
        else:
#            ignored.append(k)
            cpt +=1

    #if len(ignored):
    #    print("Ignoring : " + str(ignored) )

    return res



class DatReader(ReaderBase):
    def __init__(self):
        super(DatReader,self).__init__()
        self.commentChar= "!"
        self.readFormat = 'r'

    def ReadCleanLine(self,withError=False):
        res = super(DatReader,self).ReadCleanLine(withError=withError)
        if res is None:
            return None
        res = res.split("!")[0]

        if len(res) == 0:
            return self.ReadCleanLine(withError=withError)
        # read multiline blocks
        if res[-1] == "$":
            nline = self.PeekLine()
            if nline[0] != ".":
                return res[:-1] + " " + self.ReadCleanLine()
        return res


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

        meta = ProblemData.ProblemData()

        filetointernalid = {}
        filetointernalidElement = {}
        originalsids = []
        xs =[]
        ys = []
        zs = []
        def GetInternalNumberFromOriginalid(theList):
            try:
                return [filetointernalid[x] for x in  theList ]
            except KeyError as e:
                print("ERROR!!! Node with number " + str(e) + " not in file")
                raise

        self.FRAME  = {}
        self.AEL  = []

        l = self.ReadCleanLine()
        while True:
            if l == "RETURN":
                break
            if l is None :
                break

            ldata = LineToDic(l)

            if ldata["KEYWORD"] in KeywordToIgnore:
              self.PrintVerbose("Ignoring: " + ldata["KEYWORD"] )
              l = DischardTillNextSection(self.ReadCleanLine )
              continue

            if ldata["KEYWORD"] == "NOE":
                self.PrintVerbose('Reading Nodes')
                cpt = 0
                while(True):
                    l = self.ReadCleanLine()
                    if l is None:
                        break
                    if l[0] == ".":
                        break
                    fields = l.split()
                    if fields[0] == "I":
                        lcpt = 0
                        while lcpt < len(fields):
                            if fields[lcpt] == "I":
                                oid = int(fields[lcpt+1])
                                lcpt += 2
                            elif fields[lcpt] == "X":
                                x = float(fields[lcpt+1])
                                lcpt += 2
                            elif fields[lcpt] == "Y":
                                y = float(fields[lcpt+1])
                                lcpt += 2
                            elif fields[lcpt] == "Z":
                                z = float(fields[lcpt+1])
                                lcpt += 2
                            else:
                                raise
                    else:
                        oid = int(fields[0])
                        x = float(fields[1])
                        y = float(fields[2])
                        z = float(fields[3])
                    originalsids.append(oid)

                    xs.append(x)
                    ys.append(y)
                    zs.append(z)
                    filetointernalid[oid] = cpt
                    cpt +=1
                continue

            if ldata["KEYWORD"] == "CLM":
                self.PrintVerbose('Reading CLM')
                cpt = 0
                while(True):
                    l = self.ReadCleanLine()
                    if l is None: break
                    if l[0] == ".":
                        break
                    data = {"FIX":False,"NOEUD":False,"I":0,"C":np.zeros(3,dtype=int),"COMP":0,"NC":0.0,"V":0.0}

                    data = LineToDic(l,data)
                    if data["NOEUD"]:
                        name = "FIX_" + "".join([str(x) for x in data["C"]])
                        res.nodesTags.CreateTag(name,False).AddToTag(filetointernalid[ data["I"]])
                    elif data["COMP"] > 0 and data["V"] != 0 :
                        name = "Force"+str(data["COMP"])
                        res.nodesTags.CreateTag(name,False).AddToTag(filetointernalid[data["I"]])
                continue

            if ldata["KEYWORD"] == "SEL":
                l = l[1:] # hack to make the "while l[0] != ..." line work correctly
                fields = LineToList(l)


                while l[0] != ".":
                    cpt = 0
                    onElements = False
                    ALL = False
                    onFace = False
                    tagname = None
                    while(cpt < len(fields)):
                        if fields[cpt] == "GROUP":
                            cpt +=1
                            group = int(fields[cpt])
                            cpt +=1
                        elif fields[cpt] == "NOM":
                            cpt +=1
                            tagname = fields[cpt]
                            cpt +=1
                        elif fields[cpt] == "MAILLES":
                            onElements = True
                            cpt +=1
                        elif fields[cpt][0:6] == "NOEUD" :
                            onElements = False
                            cpt +=1
                        elif fields[cpt] == "FACES":
                            onFace = True
                            break
                        elif fields[cpt] == "TOUT":
                            ALL = True
                            cpt +=1
                        else:
                            cpt +=1
                    if onFace :
                        l = DischardTillNextSection(self.ReadCleanLine )
                        self.PrintVerbose("Ignoring Group of Faces" )
                        break
                    ids = []
                    while(True):
                        l = self.ReadCleanLine()
                        if l == None:
                            break
                        if l[0] == ".":
                            break
                        fields = LineToList(l)

                        if fields[0] == "GROUP":
                            break
                        ldata = l.split()

                        if ldata[0] == "I":
                            ids.extend(map(int,ldata[1:-1]))
                        elif ldata[0] == "UNION":
                            for name,data in res.elements.items():
                                if "Group"+str(ldata[1]) in data.tags:
                                    tag = data.tags["Group"+str(ldata[1])]
                                    ids.extend(data.originalIds[tag.GetIds()] )
                        else:
                            ids.extend(map(int,ldata))


                    if tagname is None:
                       tagname = "G"+str(group)
                    if onElements:
                        if ALL:
                            for name,data in res.elements.items():
                                data.tags.CreateTag(tagname).SetIds(np.arange(data.GetNumberOfElements()))
                        else:
                            for oidd in ids:
                                elem,idd = filetointernalidElement[oidd]
                                elem.tags.CreateTag(tagname,False).AddToTag(idd)
                                elem.tags.CreateTag("Group"+str(group),False).AddToTag(idd)
                    else:
                        if ALL:
                            res.nodesTags.CreateTag(tagname).SetIds(np.arange(len(filetointernalid) ))
                            res.nodesTags.CreateTag("Group"+str(group)).SetIds(np.arange(len(filetointernalid)))
                        else:
                            oids = [ filetointernalid[i] for i in ids]
                            res.nodesTags.CreateTag(tagname).SetIds(oids)
                            res.nodesTags.CreateTag("Group"+str(group)).SetIds(oids)
                #".SEL GROUP 1 MAILLES NOM "tagname-1_CORPS_146""
                    if l == None:
                        break
                    if l[0] == ".":
                        break

                continue

            if ldata["KEYWORD"] == "MCE":
                fields = l.split()
                if len(fields) == 1:
                    l = self.ReadCleanLine()
                    fields = l.split()

                fcpt = 0
                oid = None
                ids = None
                while(fcpt < len(fields)):
                    if fields[fcpt] == "I":
                        fcpt += 1
                        oid = int(fields[fcpt])
                        fcpt += 1
                        etype = fields[fcpt]
                        fcpt += 1
                    elif fields[fcpt] == "N":
                        fcpt += 1
                        ids = list(map(int,fields[fcpt:]))
                    else:
                        fcpt += 1
                elements = res.GetElementsOfType( EN.Bar_2)
                conn = GetInternalNumberFromOriginalid(ids[:2])
                cid = elements.AddNewElement(conn,oid)
                filetointernalidElement[oid] = (elements,cid-1)
                l = self.ReadCleanLine()

                continue

            if ldata["KEYWORD"] == "FRAME":
                self.PrintVerbose('Reading frame')
                cpt = 0
                #l = self.ReadCleanLine()
                l = self.ReadCleanLine()
                while(True): # for every I entry

                    if l == None:
                        break
                    if l[0] == ".":
                        break

                    # default values
                    __origin = np.array([0.,0.,0.])
                    __V1 = np.array([1.,0.,0.])
                    __V2 = np.array([0.,1.,0.])
                    __V3 = np.array([0.,0.,1.])
                    __TYPE = "CARTESIAN"
                    oid = -1

                    while(True):  # for every line of each I

                        if l == None:
                            break
                        if l[0] == ".":
                            break

                        fields = l.split()
                        fcpt = 0

                        while(fcpt < len(fields)): # for every field in
                            if fields[fcpt] == "I":
                                fcpt += 1
                                if oid != -1: break
                                oid = int(fields[fcpt])
                                fcpt += 1
                                continue
                            if fields[fcpt] == "TYPE":
                                fcpt += 1
                                __TYPE = (fields[fcpt])
                                fcpt += 1
                                continue
                            if fields[fcpt] == "ORIGIN":
                                fcpt += 1
                                __origin = np.array(list(map(float,fields[fcpt:fcpt+3])))
                                fcpt += 3
                                continue
                            if fields[fcpt] == "V1":
                                fcpt += 1
                                __V1 = np.array(list(map(float,fields[fcpt:fcpt+3])))
                                fcpt += 3
                                continue
                            if fields[fcpt] == "V2":
                                fcpt += 1
                                __V2 = np.array(list(map(float,fields[fcpt:fcpt+3])))
                                fcpt += 3
                                continue
                            if fields[fcpt] == "V3":
                                fcpt += 1
                                __V3 = np.array(list(map(float,fields[fcpt:fcpt+3])))
                                fcpt += 3
                                continue
                            print(fcpt)
                            print(fields)
                            raise
                        else:
                            l = self.ReadCleanLine()
                            if l == None:
                                break
                            if l[0] == ".":
                                break

                            continue

                        break

                    #print(fields)
                    self.FRAME[oid] = {"T":__TYPE, "O":__origin, "V1":__V1, "V2":__V2 }
                    #print("self.FRAME " ,self.FRAME)
                    if l == None:
                        break
                    if l[0] == ".":
                        break

                #print(self.FRAME)
                continue

            if ldata["KEYWORD"] == "AEL":
                self.PrintVerbose('Reading Elements')
                cpt = 0
                while(True):
                    l = self.ReadCleanLine()
                    if l == None:
                        break
                    if l[0] == ".":
                        break
                    fields = l.split()
                    fcpt = 0
                    __I = -1
                    __FRAME = -1
                    __MAT = -1
                    __GROUP = -1

                    while(fcpt < len(fields)):
                        if fields[fcpt] == "I":
                            fcpt += 1
                            __I = int(fields[fcpt])
                            fcpt += 1
                            continue
                        if fields[fcpt] == "FRAME":
                            fcpt += 1
                            __FRAME = int(fields[fcpt])
                            fcpt += 1
                            continue
                        if fields[fcpt] == "MAT":
                            fcpt += 1
                            __MAT = int(fields[fcpt])
                            fcpt += 1
                            continue
                        if fields[fcpt] == "GROUP":
                            fcpt += 1
                            __GROUP = fields[fcpt]
                            fcpt += 1
                            continue
                        if fields[fcpt] == "DEGRE":
                            fcpt += 1
                            __DEGRE = fields[fcpt]
                            fcpt += 1
                            continue


                        print(fields)
                        raise
                    self.AEL.append({'I':__I,"G":__GROUP,"M":__MAT,"F":__FRAME})
                continue

            if ldata["KEYWORD"] == "MAI":
                self.PrintVerbose('Reading Elements')
                #"I 1 N 55175 65855 57080 0 58679"

                cpt = 0
                while(True):
                    l = self.ReadCleanLine()
                    if l == None:
                        break
                    if l[0] == ".":
                        break

                    fields = l.split()
                    fcpt = 0
                    while(fcpt < len(fields)):
                        if fields[fcpt] == "I":
                            fcpt += 1
                            oid = int(fields[fcpt])
                            fcpt += 1
                        if fields[fcpt] == "ATT":
                            fcpt += 1
                            __data = (fields[fcpt])
                            fcpt += 1
                        if fields[fcpt] == "ENOM":
                            fcpt += 1
                            __data = (fields[fcpt])
                            fcpt += 1
                        if fields[fcpt] == "N":
                            fcpt += 1



                            p2 = []
                            while (fcpt < len(fields)):
                                p = []
                                while (fcpt < len(fields)):
                                    if int(fields[fcpt]) == 0:
                                        fcpt += 1
                                        break
                                    p.append(int(fields[fcpt]) )
                                    fcpt += 1
                                p2.append(p)
                            #self.PrintVerbose(p2)
                            if len(p2) == 2:
                                if len(p2[0]) == 3 and len(p2[1]) == 1 :
                                    elements = res.GetElementsOfType( EN.Tetrahedron_4)
                                elif len(p2[0]) == 4 and len(p2[1]) == 4 :
                                    elements = res.GetElementsOfType( EN.Hexaedron_8)
                                elif len(p2[0]) == 3 and len(p2[1]) == 3 :
                                    elements = res.GetElementsOfType( EN.Wedge_6)
                                elif len(p2[0]) == 3 and len(p2[1]) == 3 :
                                    elements = res.GetElementsOfType( EN.Wedge_6)
                                else:
                                    raise(Exception("type of element no coded" + str(p2) ))

                            elif len(p2) == 1:
                                if len(p2[0]) == 3  :
                                    elements = res.GetElementsOfType( EN.Triangle_3)
                                elif len(p2[0]) == 4 :
                                    elements = res.GetElementsOfType( EN.Quadrangle_4)
                                else:
                                    raise(Exception("type of element no coded" + str(p2) ))
                            else:
                                raise(Exception("type of element no coded" + str(p2) ))

                            flist = [item for sublist in p2 for item in sublist]
                            conn = GetInternalNumberFromOriginalid(flist)

                            cid = elements.AddNewElement(conn,oid)

                            filetointernalidElement[oid] = (elements,cid-1)

                continue
            self.SetVerboseLevel(5)
            self.PrintVerbose("--------------- Error reading this line ------------------")
            self.PrintVerbose(l)
            self.PrintVerbose("----------------------------------------------------------")
            raise

        res.SetNodes(np.array([xs,ys,zs],dtype=PBasicFloatType).T)
        res.originalIDNodes = np.array(originalsids,dtype=PBasicIndexType)
        res.PrepareForOutput()

        self.output = res

        self.filetointernalidElement = filetointernalidElement
        self.filetointernalid = filetointernalid

        #print(self.FRAME)
        #print(self.AEL)

        if len(self.FRAME) >0 and len(self.AEL)> 0:
            V1 = np.zeros( (res.GetNumberOfElements(),3))
            V1[:,0] = 1
            V2 = np.zeros( (res.GetNumberOfElements(),3))
            V2[:,1] = 1
            MAT = np.zeros( (res.GetNumberOfElements(),1))-1

            for ael  in self.AEL:
                if ael["I"] != -1:
                    oid = ael["I"]
                    oids = [ el.globaloffset+cpt  for el,cpt in [filetointernalidElement[oid]] ]

                elif ael["G"] != -1:
                    oids = res.GetElementsInTag("Group"+str(ael["G"]))
                else:
                    raise

                if ael["F"] != -1:
                    frame = self.FRAME[ael["F"]]

                    (elements,cpt) = filetointernalidElement[oid]
                    if frame["T"] != "CARTESIAN":
                        raise
                    if np.any(frame["O"]-[0,0,0]):
                        raise
                        #self.FRAME[oid] = {"T":__TYPE, "O":__origin, "V1":__V1, "V2":__V2 }
                    #print(frame)
                    V1[oids,:] = frame["V1"]
                    V2[oids,:] = frame["V2"]
                elif ael["M"] != -1:
                    MAT[oids] = ael["M"]

            res.elemFields["V1"] = V1
            res.elemFields["V2"] = V2
            res.elemFields["Mat"] = MAT
            #print(V1)
            #print(V2)


        return res


from BasicTools.IO.IOFactory import RegisterReaderClass
RegisterReaderClass(".dat",DatReader)
RegisterReaderClass(".datt",DatReader)



def CheckIntegrity():

    data = u""".INIT &
.ASEF &
MODE IMPRES 0 LECT 132 MOUCHARD 1 ECHO 1
.NOE
         1     0     0    0
         2     1.000000000E+01     0.000000000E+01     0.000000000E+00
         3     0.000000000E+01     1.000000000E+01     0.000000000E+00
         4     0.000000000E+01     0.000000000E+01     1.000000000E+00
.MAI
    I 1 N 1 2 3 0 4
.SEL GROUP 1 MAILLES NOM "SYS-1_CORPS_146"
 I 1
.SEL GROUP 2 MAILLES TOUT NOM "ALL_ELEMENTS"
.SEL GROUP 3 NOEUD TOUT NOM "ALL_NODES"
.CLM
!*********** Fixed Supports ***********
  FIX NOEUD I 1 C 1 2 3
  FIX NOEUD I 3 C 1 2 3
.FRAME
     I 1 TYPE CARTESIAN ORIGIN 0.000000 0.000000 0.000000
      V1 0.886548 -0.462637 -0.000000
      V2 0.000000 -0.000000 1.000000
     I 2 TYPE CARTESIAN ORIGIN 0.000000 0.000000 0.000000
      V1 0.959641 -0.281229 -0.000000
      V2 0.000000 -0.000000 1.000000
.AEL
     I 1 FRAME 1
"""

    res = DatReader().Read(string=data)

    from BasicTools.Helpers.Tests import TestTempDir
    tempdir = TestTempDir.GetTempPath()
    f =open(tempdir+"SamcefReader_test_File.dat","w")
    f.write(data)
    f.close()
    res = DatReader().Read(fileName=tempdir+"SamcefReader_test_File.dat")
    print(res)
    return 'ok'

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
