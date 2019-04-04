# -*- coding: utf-8 -*-
""" ASC file reader
"""

import numpy as np

__author__ = "Felipe Bordeu"

from BasicTools.Helpers.BaseOutputObject import BaseOutputObject

import BasicTools.Containers.ElementNames as EN
import BasicTools.Containers.UnstructuredMesh as UM
from BasicTools.IO.ReaderBase import ReaderBase



AscNumber = {}

AscNumber['2006'] = EN.Triangle_6
AscNumber['3010'] = EN.Tetrahedron_10
AscNumber['1002'] = EN.Bar_2

def ReadAsc(fileName=None,string=None,out=None,**kwargs):
    reader = AscReader()
    reader.SetFileName(fileName)
    reader.SetStringToRead(string)
    reader.Read(fileName=fileName, string=string,out=out,**kwargs)
    return reader.output


class AscReader(ReaderBase):
    def __init__(self):
        super(AscReader,self).__init__()
        self.commentChar= "%"
        self.readFormat = 'r'

    def Read(self,fileName=None,string=None, out=None):

        if fileName is not None:
          self.SetFileName(fileName)

        if string is not None:
          self.SetStringToRead(string)

        self.StartReading()

        if out is None:
            res = UM.UnstructuredMesh()
        else:
            res = out


        import shlex
#
#        if fileName is not None:
#            string = open(fileName, 'r')
#        elif string is not None:
#            from io import StringIO
#            string = StringIO(string)
#
#        res = UM.UnstructuredMesh()

        filetointernalid = {}
        #filetointernalidElem =  {}
        BaseOutputObject().Print("Reading file : {} ".format(fileName))

        while True:
            #line = string.readline()
            l = self.ReadCleanLine()
            if not l:
                break

            if l.find("BEGIN_NODES")>-1 :

                nbNodes = int(l.split()[1])
                BaseOutputObject().Print("Reading "+str(nbNodes)+ " Nodes")
                dim = int(l.split()[2])
                res.nodes = np.empty((nbNodes,dim))
                res.originalIDNodes= np.empty((nbNodes,))
                cpt =0;
                while(True):
                    l = self.ReadCleanLine()
                    if l.find("END_NODES") > -1:
                        break
                    s = l.split()
                    #print(s)
                    #print(res.originalIDNodes)
                    oid = int(s[0])
                    filetointernalid[oid] = cpt
                    res.originalIDNodes[cpt] = int(s[0])
                    res.nodes[cpt,:] = list(map(float,s[6:]))
                    cpt +=1
                continue

            if l.find("BEGIN_ELEMENTS")>-1 :

                nbElements = int(l.split()[1])
                BaseOutputObject().Print("Reading "+str(nbElements)+ " Elements")
                #res.nodes = np.empty((nbNodes,dim))
                #res.originalIDNodes= np.empty((nbNodes,))
                cpt =0;
                while(True):
                    l = self.ReadCleanLine()
                    if l.find("END_ELEMENTS") > -1:
                        if nbElements != cpt:# pragma: no cover
                            print("File problem!! number of elements read not equal to the total number of elemetns")
                            print(nbElements)
                            print(cpt)
                        break
                    s = l.split()

                    nametype = AscNumber[s[1]];

                    conn = [filetointernalid[x] for x in  map(int,s[5:]) ]

                    # for some types we need permutation
                    if nametype == EN.Triangle_6:
                        conn = [ conn[per] for per in [0, 2, 4, 1, 3, 5] ]
                    elif nametype == EN.Tetrahedron_10:
                        conn = [ conn[per] for per in [0, 2, 4, 9, 1,3 ,5,6,7,8] ]

                    elements = res.GetElementsOfType(nametype)
                    elements.Reserve(nbElements)

                    oid = int(s[0])

                    elements.AddNewElement(conn,oid)
                    cpt +=1
                for etype, data in res.elements.items():
                    data.tighten()
                continue

            if l.find("BEGIN_GROUPS")>-1 :

                nbgroups = int(l.split()[1])
                BaseOutputObject().Print("Reading "+str(nbgroups)+ " Groups")
                #res.nodes = np.empty((nbNodes,dim))
                #res.originalIDNodes= np.empty((nbNodes,))
                cpt =0;
                while(True):
                    l = self.ReadCleanLine()
                    if l.find("END_GROUPS") > -1:
                        if(nbgroups != (cpt)):# pragma: no cover
                            print("File problem!! number of groups read not equal to the total number of groups")
                            print(nbgroups)
                            print(cpt)
                        break
                    s = shlex.split(l)
                    tagname = s[1]
                    BaseOutputObject().Print("Reading Group " + tagname)

                    if s[2] == '1' :
                        #node group
                        tag = res.GetNodalTag(tagname)
                        tag.SetIds(np.array( [filetointernalid[x] for x in  map(int,s[7:]) ] ,dtype=np.int))
                        #tag.ids = np.zeros(len(s[7:]),dtype=np.int)
                        #cpt =0
                        #for i in s[7:]:
                        #    tag.ids[cpt] = filetointernalid[int(i)]
                        #    cpt +=1
                        #
                    else:
                        #element group

                        for x in range(7,len(s)) :
                            Oid = int(s[x])
                            #print(Oid)
                            res.AddElementToTagUsingOriginalId(Oid,tagname)
                    cpt +=1
                continue
            BaseOutputObject().PrintVerbose("Ignoring line : '" + str(l) + "'")
        self.EndReading()
        res.PrepareForOutput()
        self.output = res
        return res

from BasicTools.IO.IOFactory import RegisterReaderClass
RegisterReaderClass(".asc",AscReader)

def CheckIntegrity():

    __checkintegritydata = u"""
BEGIN_NODES 6 3
1 0 0 0 0 0 295.673175860532 28.0704731415872 346.109138100075
2 0 0 0 0 0 295.105225 28.260575 345.628395
3 0 0 0 0 0 295.180501114015 25.8084581250318 344.876373186428
4 0 0 0 0 0 295.3886425 28.1693925 345.8617875
5 0 0 0 0 0 295.231426751792 27.0671220355891 345.153077365196
6 0 0 0 0 0 295.629817604236 26.9160040797623 345.45435766729
END_NODES
BEGIN_ELEMENTS 2
21 2006 0 0 0  1 2 3 4 5 6
22 3010 0 0 0  1 2 3 1 2 3 1 2 3 1
END_ELEMENTS
BEGIN_GROUPS 2
1 PointGroug 1 0 "PART_ID 2"  ""  "PART built in Visual-Environment" 1 2 3
2 M2D 2 0 "PART_ID 2"  ""  "PART built in Visual-Environment" 21
END_GROUPS

garbage
"""
    #check from string
    res = ReadAsc(string=__checkintegritydata)
    print(res)
    if res.GetElementsOfType(EN.Triangle_6).originalIds[0] != 21:# pragma: no cover
        raise

    #check from file
    from BasicTools.Helpers.Tests import TestTempDir
    newFileName = TestTempDir().GetTempPath()+"AscFile"
    open(newFileName,'w').write(__checkintegritydata)
    res = ReadAsc(fileName = newFileName)
    print(res)
    if res.GetElementsOfType(EN.Triangle_6).originalIds[0] != 21:# pragma: no cover
        raise

    return 'ok'

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
