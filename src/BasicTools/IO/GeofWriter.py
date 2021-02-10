# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

"""Geof file writer (Zset mesh files)

"""

import numpy as np

import BasicTools.Containers.ElementNames as EN

from BasicTools.IO.WriterBase import WriterBase as WriterBase

GeofName = {}
GeofSetName= {}


from BasicTools.IO.GeofReader import PermutationZSetToBasicTools
PermutationBasicToolsToZSet = {key:np.argsort(value) for key, value in PermutationZSetToBasicTools.items() }

#0d
GeofName[EN.Point_1] = "l2d1"

#1d
GeofName[EN.Bar_2] = "l2d2"
GeofSetName[EN.Bar_2] = "line"
GeofSetName[EN.Bar_3] = "quad"


#2d
GeofName[EN.Triangle_3] = "c2d3"
GeofSetName[EN.Triangle_3] = "t3"


GeofName[EN.Triangle_6] = "c2d6"
GeofSetName[EN.Triangle_6] = "t6"
GeofSetName[EN.Quadrangle_4] = "q4"

GeofName[EN.Quadrangle_8] = "c3d8"
GeofSetName[EN.Quadrangle_8] = "q8"

#3d
GeofName[EN.Tetrahedron_4] = "c3d4"
GeofName[EN.Tetrahedron_10] = "c3d10"

GeofName[EN.Quadrangle_4] = "c2d4"
GeofName[EN.Hexaedron_8] = "c3d8"
GeofName[EN.Hexaedron_20] = "c3d20"

GeofName[EN.Wedge_6] = "c3d6"

def WriteMeshToGeof(filename,mesh, useOriginalId=False,lowerDimElementsAsElsets=False):
    OW = GeofWriter()
    OW.Open(filename)
    OW.SetWriteLowerDimElementsAsSets(lowerDimElementsAsElsets)
    OW.Write(mesh,useOriginalId = useOriginalId)
    OW.Close()

class GeofWriter(WriterBase):
    def __init__(self):
        super(GeofWriter,self).__init__()
        self.lowerDimElementsAsElsets= False
        self.canHandleBinaryChange = False

    def __str__(self):
        res  = 'GeofWriter : \n'
        res += '   FileName : '+ str(self.fileName) +'\n'
        return res

    def SetWriteLowerDimElementsAsSets(self,val):
        if val is not None:
            self.lowerDimElementsAsElsets= bool(val)

    def Write(self,meshObject,useOriginalId=False,lowerDimElementsAsElsets=None,PointFieldsNames=None,PointFields=None,CellFieldsNames=None,CellFields=None):

        if PointFieldsNames is not None or \
           PointFields      is not None or \
           CellFieldsNames  is not None or \
           CellFields       is not None:
               print("warning GeofWriter only canr write the mesh, fields are ignored")

        self.SetWriteLowerDimElementsAsSets(lowerDimElementsAsElsets)

        meshObject.PrepareForOutput()

        self.filePointer.write("% Written by BasicTools package\n")

        self.filePointer.write("***geometry\n");
        self.filePointer.write("  **node\n");
        numberofpoints = meshObject.GetNumberOfNodes()
        self.filePointer.write("{} {} \n".format(numberofpoints,meshObject.GetDimensionality()) )
        #
        posn = meshObject.GetPosOfNodes()
        if useOriginalId:
           for n in range(numberofpoints):
               self.filePointer.write("{} ".format(int(meshObject.originalIDNodes[n])))
               #np.savetxt(self.filePointer, posn[n,:] )
               posn[n,:].tofile(self.filePointer, sep=" ")
               self.filePointer.write("\n")
        else:
           for n in range(numberofpoints):
               self.filePointer.write("{} ".format(n+1) )
               #np.savetxt(self.filePointer, posn[np.newaxis,n,:] )
               posn[np.newaxis,n,:].tofile(self.filePointer, sep=" ")
               self.filePointer.write("\n")
        #
        nbElements = 0
        maxDimensionalityOfelements = 0
        for ntype,elems in meshObject.elements.items():
            maxDimensionalityOfelements = max(EN.dimension[ntype],maxDimensionalityOfelements)

        for ntype,elems in meshObject.elements.items():
            if EN.dimension[ntype] == maxDimensionalityOfelements or self.lowerDimElementsAsElsets == True:
                nbElements += elems.GetNumberOfElements()


        self.filePointer.write("  **element\n")
        self.filePointer.write("{}\n".format(nbElements))


        cpt =0;
        for ntype, data in meshObject.elements.items():
            elemtype = GeofName[ntype]
            if EN.dimension[ntype] != maxDimensionalityOfelements and self.lowerDimElementsAsElsets == False:
                continue
            #npe = data.GetNumberOfNodesPerElement()
            #if elemtype!="c2d3":
            for i in range(data.GetNumberOfElements() ):
                conn = data.connectivity[i,:].ravel()

                if elemtype in PermutationBasicToolsToZSet:
                  conn = [conn[x] for x in PermutationBasicToolsToZSet[elemtype]]
                if useOriginalId:
                    self.filePointer.write("{} {} ".format(data.originalIds[i],elemtype) )
                    self.filePointer.write(" ".join([str(int(meshObject.originalIDNodes[x])) for x in conn]))
                else:
                    self.filePointer.write("{} {} ".format(cpt+1,elemtype) )
                    self.filePointer.write(" ".join([str(x+1) for x in conn]))
                cpt += 1;
                self.filePointer.write("\n")

        self.filePointer.write(" ***group \n")

        for tag in meshObject.nodesTags:
            if len(tag) == 0:
                continue
            self.filePointer.write("  **nset {} \n".format(tag.name))
            data = np.zeros((meshObject.GetNumberOfNodes(),1),dtype=np.int)
            if useOriginalId:
                self.filePointer.write(" ".join([str(int(meshObject.originalIDNodes[x])) for x in tag.GetIds()]))
            else:
                self.filePointer.write(" ".join([str(x+1) for x in tag.GetIds()]))
            self.filePointer.write("\n")

        meshObject.PrepareForOutput();

        if self.lowerDimElementsAsElsets == False :
            celtags = meshObject.GetNamesOfElemTags()
            for tagname in celtags:

                idInTag = 0
                flag = False
                for ntype,elems in meshObject.elements.items():
                    if EN.dimension[ntype] == maxDimensionalityOfelements:
                        if tagname in elems.tags:
                            flag =  True
                            idInTag += elems.tags[tagname].cpt
                self.PrintVerbose("Tag " + str(tagname) + " has "+ str(idInTag) + " elements")
                # no output if no elements in tag
                if flag == False : continue
                #empty tags
                if idInTag == 0 :  continue
                self.filePointer.write("  **elset {} \n".format(tagname))
                cpt =0

                for ntype,elems in meshObject.elements.items():
                    if EN.dimension[ntype] != maxDimensionalityOfelements:
                        continue
                    if tagname in elems.tags:
                        tag = elems.tags[tagname]
                        if tag.cpt :
                            self.filePointer.write(" ".join([str(x+1+cpt) for x in tag.GetIds()]))
                    cpt += elems.GetNumberOfElements()

                self.filePointer.write("\n")
        else:
            celtags = meshObject.GetNamesOfElemTags()
            for tagname in celtags:
                self.filePointer.write("  **elset {} \n".format(tagname))
                data = meshObject.GetElementsInTag(tagname,useOriginalId=useOriginalId)
                if useOriginalId :
                    self.filePointer.write(" ".join([str(x) for x in data]))
                else:
                    self.filePointer.write(" ".join([str(x+1) for x in data]))
                self.filePointer.write("\n")

        # Dotsets, lisets, facets
        if self.lowerDimElementsAsElsets == False:
            for dimToTreat in range(maxDimensionalityOfelements):

                celtags = meshObject.GetNamesOfElemTags()
                for tagname in celtags:

                    idInTag = 0
                    flag = False
                    for ntype,elems in meshObject.elements.items():
                        if EN.dimension[ntype] == dimToTreat:
                            if tagname in elems.tags:
                                flag =  True
                                idInTag += elems.tags[tagname].cpt
                    self.PrintVerbose("Set  " + str(tagname) + " has "+ str(idInTag) + " elements")
                    # no output if no elements in tag
                    if flag == False : continue
                    #empty tags
                    if idInTag == 0 :  continue
                    if dimToTreat == 0:# pragma: no cover
                        raise Exception("Dont know how to treat the doset, please code me")
                    #    self.filePointer.write("  **doset {} \n".format(tagname))
                    elif dimToTreat == 1:
                        self.filePointer.write("  **liset {} \n".format(tagname))
                    elif dimToTreat == 2:
                        self.filePointer.write("  **faset {} \n".format(tagname))


                    for ntype,elems in meshObject.elements.items():
                        if EN.dimension[ntype] != dimToTreat:
                            continue
                        if tagname in elems.tags:
                            tag = elems.tags[tagname]

                            name = GeofSetName[ntype];

                            ids = tag.GetIds()
                            for e in range(len(tag)):
                                conn = elems.connectivity[ids[e],:]
                                if name in PermutationBasicToolsToZSet:
                                    conn = [conn[x] for x in PermutationBasicToolsToZSet[name]]
                                self.filePointer.write(" {} ".format(name))
                                self.filePointer.write(" ".join([str(x+1) for x in conn ]))
                                self.filePointer.write(" \n")

        self.filePointer.write("***return \n")

from BasicTools.IO.IOFactory import RegisterWriterClass
RegisterWriterClass(".geof",GeofWriter)

def CheckIntegrity():
    import BasicTools.Containers.UnstructuredMesh as UM

    from BasicTools.Helpers.Tests import TestTempDir

    tempdir = TestTempDir.GetTempPath()
    print(tempdir)

    mymesh = UM.UnstructuredMesh()
    mymesh.nodes = np.array([[0.00000000001,0,0],
                             [1,0,0],
                             [0,1,0],
                             [1,1,0],
                             [0.5,0,0.1],
                             [0,0.5,0.1],
                             [0.5,0.5,0.1],
    ],dtype=np.float)
    mymesh.originalIDNodes = np.array([1, 3, 4, 5, 6, 7, 8],dtype=np.int)

    mymesh.nodesTags.CreateTag("coucou").AddToTag(0)

    tet = mymesh.GetElementsOfType(EN.Tetrahedron_4)
    tet.AddNewElement([0,1,2,3],0)
    tet.tags.CreateTag("TheOnlyTet").AddToTag(0)

    tris = mymesh.GetElementsOfType(EN.Triangle_3)
    tris.AddNewElement([0,1,2],0)
    tris.AddNewElement([2,1,3],3)
    tris.originalIds = np.array([3, 5],dtype=np.int)

    bars = mymesh.GetElementsOfType(EN.Bar_2)
    bars.AddNewElement([0,1],0)
    bars.AddNewElement([1,3],1)
    bars.tags.CreateTag("firstBar").AddToTag(0)

    bars3 = mymesh.GetElementsOfType(EN.Triangle_6)
    bars3.AddNewElement([0,1,2,4,6,5],0)
    bars3.tags.CreateTag("Tri6").AddToTag(0)

    #point = mymesh.GetElementsOfType(EN.Point_1)
    #point.AddNewElement([0],0)
    #point.tags.CreateTag("onlyPoint").AddToTag(0)


    mymesh.AddElementToTagUsingOriginalId(3,"Tag1")
    mymesh.AddElementToTagUsingOriginalId(5,"Tag3")

    OW = GeofWriter()
    print(OW)
    OW.Open(tempdir+"Test_GeoWriter.geof")
    OW.Write(mymesh, useOriginalId=True)
    OW.Close()

    WriteMeshToGeof(tempdir+"Test_GeoWriter_II.geof", mymesh,useOriginalId=False )
    WriteMeshToGeof(tempdir+"Test_GeoWriter_III.geof", mymesh,useOriginalId=False,lowerDimElementsAsElsets=True )
    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
