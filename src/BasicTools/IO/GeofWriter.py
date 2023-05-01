# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

"""Geof file writer (Zset mesh files)
"""
from typing import Optional

import numpy as np

import BasicTools.Containers.ElementNames as EN

from BasicTools.IO.WriterBase import WriterBase
from BasicTools.NumpyDefs import  PBasicFloatType, PBasicIndexType
from BasicTools.Containers.UnstructuredMesh import UnstructuredMesh
from BasicTools.IO.GeofReader import PermutationZSetToBasicTools
PermutationBasicToolsToZSet = {key:np.argsort(value) for key, value in PermutationZSetToBasicTools.items() }

GeofName = {}
GeofSetName= {}

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

def WriteMeshToGeof(fileName :str , mesh:UnstructuredMesh, useOriginalId:Optional[bool]=False,lowerDimElementsAsElsets:bool=False)-> None:
    """Function API for writing mesh in the geof format file.
        A file is created using the path and name of fileName

    Parameters
    ----------
    fileName : str
        name with path to the file to be created (relative or absolute)
    mesh : UnstructuredMesh
        the mesh to be exported
    useOriginalId : Optional[bool], optional
        use the Original Id for the number of nodes and elements
        (the user is responsible of the consistency of this data), by default False
    lowerDimElementsAsElsets : Optional[bool], optional
        refert to GeofWriter.SetWriteLowerDimElementsAsSets for the documentations
    """

    geofWriter = GeofWriter()
    geofWriter.Open(fileName)
    geofWriter.SetWriteLowerDimElementsAsSets(lowerDimElementsAsElsets)
    geofWriter.Write(mesh,useOriginalId = useOriginalId)
    geofWriter.Close()

class GeofWriter(WriterBase):
    """Class to write Unstructured mesh on disk in the geof format file
    """
    def __init__(self) -> None:
        super(GeofWriter,self).__init__()
        self.lowerDimElementsAsElsets= False
        self.canHandleBinaryChange = False

    def __str__(self):
        res  = 'GeofWriter : \n'
        res += '   FileName : '+ str(self.fileName) +'\n'
        return res

    def SetWriteLowerDimElementsAsSets(self,val:bool) -> None:
        """ Set the type of export for element of lower dimensionality
            (1D and 2D element on a 3D mesh)

        Parameters
        ----------
        val : bool
            if True all the elements of lower dimensionality are exported in the mesh
            if False all the element of lower dimensionality are exported as bsets,
            by default False
        """
        if val is not None:
            self.lowerDimElementsAsElsets= bool(val)

    def Write(self,meshObject: UnstructuredMesh, useOriginalId:Optional[bool]=False, lowerDimElementsAsElsets:bool=False, PointFieldsNames=None, PointFields=None, CellFieldsNames=None, CellFields=None):
        """Write mesh to file in Geof format

        Parameters
        ----------
        meshObject : UnstructuredMesh
            the mesh to be written
        useOriginalId : Optional[bool], optional
            use the Original Id for the number of nodes and elements
            (the user is responsible of the consistency of this data),
            by default False
        lowerDimElementsAsElsets : Optional[bool], optional
            refert to GeofWriter.SetWriteLowerDimElementsAsSets for the documentations
        PointFieldsNames : None
            Not Used, by default None
        PointFields : None
            Not Used, by default None
        CellFieldsNames : None
            Not Used, by default None
        CellFields : None
            Not Used, by default None

        Raises
        ------
        Exception
            in the case 1D element are presents in the mesh (doset)
        """
        if PointFieldsNames is not None or PointFields      is not None or \
            CellFieldsNames  is not None or CellFields       is not None:
            print("warning GeofWriter only can write the mesh, fields are ignored")

        self.SetWriteLowerDimElementsAsSets(lowerDimElementsAsElsets)

        meshObject.PrepareForOutput()

        self.filePointer.write("% Written by BasicTools package\n")

        self.filePointer.write("***geometry\n")
        self.filePointer.write("  **node\n")
        numberOfPoints = meshObject.GetNumberOfNodes()
        self.filePointer.write(f"{numberOfPoints} {meshObject.GetDimensionality()} \n" )
        #
        posn = meshObject.GetPosOfNodes()
        if useOriginalId:
            for n in range(numberOfPoints):
                self.filePointer.write(f"{int(meshObject.originalIDNodes[n])} ")
                posn[n,:].tofile(self.filePointer, sep=" ")
                self.filePointer.write("\n")
        else:
            for n in range(numberOfPoints):
                self.filePointer.write(f"{n+1} ")
                posn[np.newaxis,n,:].tofile(self.filePointer, sep=" ")
                self.filePointer.write("\n")
        #
        nbElements = 0
        maxDimensionalityOfElements = 0
        for ntype,elems in meshObject.elements.items():
            maxDimensionalityOfElements = max(EN.dimension[ntype],maxDimensionalityOfElements)

        for ntype,elems in meshObject.elements.items():
            if EN.dimension[ntype] == maxDimensionalityOfElements or self.lowerDimElementsAsElsets is True:
                nbElements += elems.GetNumberOfElements()


        self.filePointer.write("  **element\n")
        self.filePointer.write(f"{nbElements}\n")


        cpt =0
        for ntype, data in meshObject.elements.items():
            elemtype = GeofName[ntype]
            if EN.dimension[ntype] != maxDimensionalityOfElements and self.lowerDimElementsAsElsets is False:
                continue
            #npe = data.GetNumberOfNodesPerElement()
            #if elemtype!="c2d3":
            for i in range(data.GetNumberOfElements() ):
                conn = data.connectivity[i,:].ravel()

                if elemtype in PermutationBasicToolsToZSet:
                    conn = [conn[x] for x in PermutationBasicToolsToZSet[elemtype]]
                if useOriginalId:
                    self.filePointer.write(f"{data.originalIds[i]} {elemtype} " )
                    self.filePointer.write(" ".join([str(int(meshObject.originalIDNodes[x])) for x in conn]))
                else:
                    self.filePointer.write(f"{cpt+1} {elemtype} ")
                    self.filePointer.write(" ".join([str(x+1) for x in conn]))
                cpt += 1
                self.filePointer.write("\n")

        self.filePointer.write(" ***group \n")

        for tag in meshObject.nodesTags:
            if len(tag) == 0:
                continue
            self.filePointer.write(f"  **nset {tag.name} \n")
            data = np.zeros((meshObject.GetNumberOfNodes(),1),dtype=PBasicIndexType)
            if useOriginalId:
                self.filePointer.write(" ".join([str(int(meshObject.originalIDNodes[x])) for x in tag.GetIds()]))
            else:
                self.filePointer.write(" ".join([str(x+1) for x in tag.GetIds()]))
            self.filePointer.write("\n")

        meshObject.PrepareForOutput()

        if self.lowerDimElementsAsElsets is False :
            celtags = meshObject.GetNamesOfElemTags()
            for tagname in celtags:

                idInTag = 0
                flag = False
                for ntype,elems in meshObject.elements.items():
                    if EN.dimension[ntype] == maxDimensionalityOfElements:
                        if tagname in elems.tags:
                            flag =  True
                            idInTag += elems.tags[tagname].cpt
                self.PrintVerbose(f"Tag {tagname} has {idInTag} elements")
                # no output if no elements in tag
                if flag is False :
                    continue
                #empty tags
                if idInTag == 0 :
                    continue
                self.filePointer.write(f"  **elset {tagname} \n")
                cpt =0

                for ntype,elems in meshObject.elements.items():
                    if EN.dimension[ntype] != maxDimensionalityOfElements:
                        continue
                    if tagname in elems.tags:
                        tag = elems.tags[tagname]
                        if tag.cpt :
                            if useOriginalId:
                                self.filePointer.write(" ".join([str(elems.originalIds[x]) for x in tag.GetIds()]))
                            else:
                                self.filePointer.write(" ".join([str(x+1+cpt) for x in tag.GetIds()]))
                            self.filePointer.write(" ")
                    cpt += elems.GetNumberOfElements()

                self.filePointer.write("\n")
        else:
            celtags = meshObject.GetNamesOfElemTags()
            for tagname in celtags:
                self.filePointer.write(f"  **elset {tagname} \n")
                data = meshObject.GetElementsInTag(tagname,useOriginalId=useOriginalId)
                if useOriginalId :
                    self.filePointer.write(" ".join([str(x) for x in data]))
                    self.filePointer.write(" ")
                else:
                    self.filePointer.write(" ".join([str(x+1) for x in data]))
                    self.filePointer.write(" ")
                self.filePointer.write("\n")

        # Dotsets, lisets, facets
        if self.lowerDimElementsAsElsets is False:
            for dimToTreat in range(maxDimensionalityOfElements):

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
                    if flag is False :
                        continue
                    #empty tags
                    if idInTag == 0 :
                        continue
                    if dimToTreat == 0:# pragma: no cover
                        raise Exception("Don't know how to treat the doset, please code me")
                    #    self.filePointer.write("  **doset {} \n".format(tagname))
                    elif dimToTreat == 1:
                        self.filePointer.write(f"  **liset {tagname} \n")
                    elif dimToTreat == 2:
                        self.filePointer.write(f"  **faset {tagname} \n")


                    for ntype,elems in meshObject.elements.items():
                        if EN.dimension[ntype] != dimToTreat:
                            continue
                        if tagname in elems.tags:
                            tag = elems.tags[tagname]

                            name = GeofSetName[ntype]

                            ids = tag.GetIds()
                            for e in range(len(tag)):
                                conn = elems.connectivity[ids[e],:]
                                if name in PermutationBasicToolsToZSet:
                                    conn = [conn[x] for x in PermutationBasicToolsToZSet[name]]
                                self.filePointer.write(f" {name} ")
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
    ],dtype=PBasicFloatType)
    mymesh.originalIDNodes = np.array([1, 3, 4, 5, 6, 7, 8],dtype=PBasicIndexType)

    mymesh.nodesTags.CreateTag("coucou").AddToTag(0)
    mymesh.nodesTags.CreateTag("empty tag")

    tet = mymesh.GetElementsOfType(EN.Tetrahedron_4)
    tet.AddNewElement([0,1,2,3],0)
    tet.tags.CreateTag("TheOnlyTet").AddToTag(0)

    triangles = mymesh.GetElementsOfType(EN.Triangle_3)
    triangles.AddNewElement([0,1,2],0)
    triangles.AddNewElement([2,1,3],3)
    triangles.originalIds = np.array([3, 5],dtype=PBasicIndexType)

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

    geofWriter = GeofWriter()
    print(geofWriter)
    geofWriter.Open(tempdir+"Test_GeoWriter.geof")
    geofWriter.Write(mymesh, useOriginalId=True, PointFieldsNames=[], PointFields=[])
    geofWriter.Close()

    WriteMeshToGeof(tempdir+"Test_GeoWriter_II.geof", mymesh,useOriginalId=False )
    WriteMeshToGeof(tempdir+"Test_GeoWriter_III.geof", mymesh,useOriginalId=False,lowerDimElementsAsElsets=True )
    WriteMeshToGeof(tempdir+"Test_GeoWriter_III.geof", mymesh,useOriginalId=True,lowerDimElementsAsElsets=True )

    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
