# -*- coding: utf-8 -*-
""" Gmsh file reader (gmesh mesh files)

"""
import numpy as np
__author__ = "Felipe Bordeu"

import BasicTools.Containers.ElementNames as EN
import BasicTools.Containers.UnstructuredMesh  as UM


gmshNumber = {}
gmshNumber['1'] = EN.Bar_2
gmshNumber['2'] = EN.Triangle_3
gmshNumber['3'] = EN.Quadrangle_4
gmshNumber['4'] = EN.Tetrahedron_4
gmshNumber['5'] = EN.Hexaedron_8
gmshNumber['6'] = EN.Wedge_6
gmshNumber['7'] = EN.Pyramid_5
gmshNumber['15'] = EN.Point_1


def ReadGmsh(fileName=None,string=None,out=None):

    if fileName is not None:
        string = open(fileName, 'r')
    elif string is not None:
        from io import StringIO
        string = StringIO(string)

    if out is None:
        res = UM.UnstructuredMesh()
    else:
        res = out # pragma: no cover

    filetointernalid = {}
    #filetointernalidElem =  {}
    tagsNames = []
    while(True):
        line = string.readline()
        if not line: break

        l = line.strip('\n').lstrip().rstrip()
        if len(l) == 0: continue

        if l.find("$MeshFormat")>-1 :
            while(True):
                line = string.readline()
                l = line.strip('\n').lstrip().rstrip()
                if len(l) == 0: continue
                if l.find("$EndMeshFormat") > -1:
                    break
            continue

        if l.find("$Nodes")>-1 :
            line = string.readline()
            l = line.strip('\n').lstrip().rstrip()

            nbNodes = int(l.split()[0])
            #print("Reading "+str(nbNodes)+ " Nodes")
            res.nodes = np.empty((nbNodes,3))
            res.originalIDNodes= np.empty((nbNodes,),dtype=np.int)
            cpt =0;
            while(True):
                line = string.readline()
                l = line.strip('\n').lstrip().rstrip()
                if len(l) == 0: continue
                if l.find("$EndNodes") > -1:
                    break
                s = l.split()
                #print(s)
                #print(res.originalIDNodes)
                oid = int(s[0])
                if cpt == nbNodes :
                    raise(Exception("More points than the number of point in the header (fix your file!!!)")) # pragma: no cover
                filetointernalid[oid] = cpt
                res.originalIDNodes[cpt] = int(s[0])
                res.nodes[cpt,:] = list(map(float,s[1:]))
                cpt +=1

            continue
        if l.find("$PhysicalNames")>-1 :
            line = string.readline()
            l = line.strip('\n').strip()
            numberOfTags = int(l)
            print(numberOfTags)
            if l.find("$EndPhysicalNames") > -1:
                continue

            while(numberOfTags > 0):
                numberOfTags -= 1
                line = string.readline()
                if line.find("$EndPhysicalNames") > -1:
                    break
                l = line.strip('\n').strip().split()

                dimensionality = int(l[0])
                tagNumber = l[1]
                tagName = l[2].strip('"').strip("'")
                tagsNames.append((dimensionality,tagNumber,tagName)  )

            line = string.readline()
            if line.find("$EndPhysicalNames") == -1:
                print("Specting $EndPhysicalNames Keyword in the file (corrupted file??")

            continue


        if l.find("$Elements")>-1 :
            line = string.readline()
            l = line.strip('\n').lstrip().rstrip()

            nbElements = int(l.split()[0])
            #print("Reading "+str(nbElements)+ " Elements")
            #res.nodes = np.empty((nbNodes,dim))
            #res.originalIDNodes= np.empty((nbNodes,))
            cpt =0;
            while(True):
                line = string.readline()
                l = line.strip('\n').lstrip().rstrip()
                if len(l) == 0: continue
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

                conn = [filetointernalid[x] for x in  map(int,s[(ntags+3):]) ]
                elements = res.GetElementsOfType(nametype)
                elements.AddNewElement(conn,oid)

                if ntags >=0 :
                    res.AddElementToTagUsingOriginalId(oid,"PhyTag"+str(s[3]))
                if ntags >=1 :
                    res.AddElementToTagUsingOriginalId(oid,"GeoTag"+str(s[4]))

                for n in range(2, ntags):
                    res.AddElementToTagUsingOriginalId(oid,"ExtraTag"+str(n-2))

                cpt +=1
            continue
        print("ignoring line : " + l )


    ## apply tags names
    for dim,number,newName in tagsNames:
        print("*** chaning tag form '" + str("PhyTag"+number) + "' to " + str(newName)+  " of dim " + str(  dim) )
        for name,data in res.elements.items():
            print(dim)
            print(EN.dimension[name])
            if EN.dimension[name] != dim :
                continue
            print("chaning tag form " + str("PhyTag"+number) + " to " + str(newName))
            data.tags.RenameTag("PhyTag"+number,newName,noError=True)

    res.PrepareForOutput()
    return res



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
