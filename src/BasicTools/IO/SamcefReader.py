#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
# -*- coding: utf-8 -*-
import numpy as np

import BasicTools.Containers.ElementNames as EN
from BasicTools.IO.ReaderBase import ReaderBase
from BasicTools.Helpers import ParserHelper as PH

KeywordToIgnore = ["INIT",
                   "ASEF",
                   "HYP",
                   "AEL",
                   "DES",
                   "MAT",
                   "UNIT",
                   "GEL",
                   "SAM",
                   "OPT",
                   "SAI",
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
        #if "&" in res:
        #    res = res.split("&")[1]
        if len(res) == 0:
            return self.ReadCleanLine(withError=withError)
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
                fields = LineToList(l)


                cpt = 0
                onElements = False
                ALL = False
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
                    elif fields[cpt] == "NOEUD":
                        onElements = False
                        cpt +=1
                    elif fields[cpt] == "TOUT":
                        ALL = True
                        cpt +=1
                    else:
                        cpt +=1

                ids = []
                while(True):
                    l = self.ReadCleanLine()
                    if l == None:
                        break
                    if l[0] == ".":
                        break

                    lis = l.replace("I","").replace("$","").split()
                    ids.extend(map(int,lis))

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
                        res.nodesTags.CreateTag(tagname).SetIds(ids)
                        res.nodesTags.CreateTag("Group"+str(group)).SetIds(ids)
                #".SEL GROUP 1 MAILLES NOM "tagname-1_CORPS_146""
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
                                    #tetra
                                    elements = res.GetElementsOfType( EN.Tetrahedron_4)
                                elif len(p2[0]) == 4 and len(p2[1]) == 4 :
                                    elements = res.GetElementsOfType( EN.Hexaedron_8)
                                elif len(p2[0]) == 3 and len(p2[1]) == 3 :
                                    elements = res.GetElementsOfType( EN.Wedge_6)
                                else:
                                    raise(Exception("type of element no coded" + str(p2) ))
                                flist = [item for sublist in p2 for item in sublist]
                                conn = GetInternalNumberFromOriginalid(flist)

                                cid = elements.AddNewElement(conn,oid)
                                filetointernalidElement[oid] = (elements,cid-1)
                continue
            self.PrintVerbose(l)
            raise

        res.nodes = np.array([xs,ys,zs],dtype=np.float).T
        res.originalIDNodes = np.array(originalsids,dtype=np.int)
        res.PrepareForOutput()

        self.output = res
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
