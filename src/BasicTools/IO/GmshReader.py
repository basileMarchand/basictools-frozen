# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

""" Gmsh file reader (gmesh mesh files)

"""
import numpy as np

import BasicTools.Containers.ElementNames as EN
import BasicTools.Containers.UnstructuredMesh  as UM
from BasicTools.IO.ReaderBase import ReaderBase
from BasicTools.NumpyDefs import PBasicIndexType

gmshNumber = {}
gmshNumber['1'] = EN.Bar_2
gmshNumber['2'] = EN.Triangle_3
gmshNumber['3'] = EN.Quadrangle_4
gmshNumber['4'] = EN.Tetrahedron_4
gmshNumber['5'] = EN.Hexaedron_8
gmshNumber['6'] = EN.Wedge_6
gmshNumber['7'] = EN.Pyramid_5
gmshNumber['8'] = EN.Bar_3
gmshNumber['9'] = EN.Triangle_6
gmshNumber['10'] = EN.Quadrangle_9
gmshNumber['15'] = EN.Point_1
gmshNumber['16'] = EN.Quadrangle_8
gmshNumber['17'] = EN.Hexaedron_20


def ReadGmsh(fileName=None,string=None,out=None,**kwargs):
    reader = GmshReader()
    reader.SetFileName(fileName)
    reader.SetStringToRead(string)
    reader.Read(fileName=fileName, string=string,out=out,**kwargs)
    return reader.output


class GmshReader(ReaderBase):
    def __init__(self):
        super(GmshReader,self).__init__()
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


        filetointernalid = {}
        #filetointernalidElem =  {}
        tagsNames = []
        while(True):
            l = self.ReadCleanLine()
            if not l: break

            if l.find("$MeshFormat")>-1 :
                while(True):
                    line = self.ReadCleanLine()
                    l = line.strip('\n').lstrip().rstrip()
                    if len(l) == 0: continue
                    if l.find("$EndMeshFormat") > -1:
                        break
                continue

            if l.find("$Nodes")>-1 :
                l = self.ReadCleanLine()

                nbNodes = int(l.split()[0])
                #print("Reading "+str(nbNodes)+ " Nodes")
                res.nodes = np.empty((nbNodes,3))
                res.originalIDNodes= np.empty((nbNodes,),dtype=PBasicIndexType)
                cpt =0;
                while(True):
                    line = self.ReadCleanLine()
                    l = line.strip('\n').lstrip().rstrip()
                    if len(l) == 0: continue
                    if l.find("$EndNodes") > -1:
                        break
                    s = l.split()
                    if cpt == nbNodes :
                        raise(Exception("More points than the number of point in the header (fix your file!!!)")) # pragma: no cover
                    filetointernalid[s[0]] = cpt
                    res.originalIDNodes[cpt] = int(s[0])
                    res.nodes[cpt,:] = list(map(float,s[1:]))
                    cpt +=1

                continue
            if l.find("$PhysicalNames")>-1 :
                l = self.ReadCleanLine()
                numberOfTags = int(l)
                if l.find("$EndPhysicalNames") > -1:
                    continue

                while(numberOfTags > 0):
                    numberOfTags -= 1
                    l = self.ReadCleanLine()
                    if l.find("$EndPhysicalNames") > -1:
                        break
                    l = l.split()

                    dimensionality = int(l[0])
                    tagNumber = l[1]
                    tagName = l[2].strip('"').strip("'")
                    tagsNames.append((dimensionality,tagNumber,tagName)  )

                line = self.ReadCleanLine()
                if line.find("$EndPhysicalNames") == -1:
                    print("Specting $EndPhysicalNames Keyword in the file (corrupted file??")

                continue


            if l.find("$Elements")>-1 :
                line = self.ReadCleanLine()
                l = line.strip('\n').lstrip().rstrip()

                nbElements = int(l.split()[0])
                #print("Reading "+str(nbElements)+ " Elements")
                #res.nodes = np.empty((nbNodes,dim))
                #res.originalIDNodes= np.empty((nbNodes,))
                allTags = {}
                cpt =0;
                while(True):
                    l = self.ReadCleanLine()
                    if l.find("$EndElements") > -1:
                        if nbElements != cpt:# pragma: no cover
                            print("File problem!! number of elements read not equal to the total number of elements")
                            print(nbElements)
                            print(cpt)
                        break
                    s = l.split()

                    oid = int(s[0])
                    gmshElemType = s[1]
                    nametype = gmshNumber[gmshElemType]

                    ntags = int(s[2])

                    conn = [filetointernalid[x] for x in  s[ntags+3:] ]
                    elements = res.elements.GetElementsOfType(nametype)
                    elnumber = elements.AddNewElement(conn,oid)-1
                    tags = allTags.setdefault(nametype,{})
                    if ntags >=0 :
                        tags.setdefault("PhyTag"+s[3],[]).append(elnumber)
                    if ntags >=1 :
                        tags.setdefault("GeoTag"+s[4],[]).append(elnumber)
                    for n in range(2, ntags):
                        tags.setdefault("ExtraTag"+str(n-2),[]).append(elnumber)

                    cpt +=1
                for nametype,tags in allTags.items():
                    elements = res.elements.GetElementsOfType(nametype)
                    for name,ids in tags.items():
                        elements.tags.CreateTag(name).SetIds(ids)
                del allTags
                continue
            print("ignoring line : " + l )


        ## convert to 2D if needed
        if np.linalg.norm(res.nodes[:,2]) == 0:
            from BasicTools.Containers.UnstructuredMeshModificationTools import LowerNodesDimension
            res = LowerNodesDimension(res)


        ## apply tags names

        for dim,number,newName in tagsNames:
            for name,data in res.elements.items():
                if EN.dimension[name] != dim :
                    continue
                self.PrintVerbose("changing tag form '" + str("PhyTag"+number) + "' to '" + str(newName)+"'")
                data.tags.RenameTag("PhyTag"+number,newName,noError=True)

        self.EndReading()
        res.PrepareForOutput()
        self.output = res
        return res

from BasicTools.IO.IOFactory import RegisterReaderClass
RegisterReaderClass(".msh",GmshReader)



def CheckIntegrity():

    __teststring = u"""
$MeshFormat
2.2 0 8
$EndMeshFormat
$PhysicalNames
2
0 223 "1D"
0 227 "2D"
$EndPhysicalNames
$Nodes
3
1 30 0 0
2 30 0 75
3 30 -2.5 0
$EndNodes
$Elements
2
1 15 2 223 1 2
2 15 3 227 2 3 1
$EndElements
this is a comment
"""

    res = ReadGmsh(string=__teststring)

    print("----")
    print(res.nodes)
    print(res.originalIDNodes)
    print(res.GetElementsOfType('bar2').connectivity)


    from BasicTools.Helpers.Tests import TestTempDir
    newFileName = TestTempDir().GetTempPath()+"mshFile"
    open(newFileName,'w').write(__teststring)
    res = ReadGmsh(fileName=newFileName)
    print(res)
    return 'ok'

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
