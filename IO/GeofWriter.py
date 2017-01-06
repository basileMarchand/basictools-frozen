# -*- coding: utf-8 -*-

from OTTools.FE.MeshBase import Tag as Tag
from OTTools.IO.WriterBase import WriterBase as WriterBase
import numpy as np
import OTTools.FE.ElementNames as EN

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
#3d
GeofName[EN.Tetrahedron_4] = "c3d4"
GeofName[EN.Tetrahedron_10] = "c3d10"
GeofName[EN.Quadrangle_4] = "c2d4"
GeofName[EN.Hexaedron_8] = "c3d8"

GeofName[EN.Wedge_6] = "c3d6"

def WriteMeshToGeof(filename,mesh, useOriginalId=False,lowerDimElementsAsSets=False):
    OW = GeofWriter()
    OW.Open(filename)
    OW.Write(mesh,useOriginalId = useOriginalId,lowerDimElementsAsSets=lowerDimElementsAsSets)
    OW.Close()

class GeofWriter(WriterBase):
    def __init__(self):
        super(GeofWriter,self).__init__()

    def __str__(self):
        res  = 'GeofWriter : \n'
        res += '   FileName : '+ str(self.fileName) +'\n'
        return res

    def SetFileName(self,fileName):
        self.fileName = fileName;

    def Write(self,meshObject,useOriginalId=False,lowerDimElementsAsSets=False):

        self.filePointer.write("% This file has been writen by the python routine GeofWriter of the OTTools package\n")
        self.filePointer.write("% For any question about this routine, please contact SAFRAN TECH Pole M&S Team OT\n")

        self.filePointer.write("***geometry\n");
        self.filePointer.write("  **node\n");
        numberofpoints = meshObject.GetNumberOfNodes()
        self.filePointer.write("{} {} \n".format(numberofpoints,3) )
        #
        posn = meshObject.GetPosOfNodes()
        if useOriginalId:
           for n in xrange(numberofpoints):
               self.filePointer.write("{} ".format(int(meshObject.originalIDNodes[n])))
               np.savetxt(self.filePointer, posn[np.newaxis,n,:] )
        else:
           for n in xrange(numberofpoints):
               self.filePointer.write("{} ".format(n+1) )
               np.savetxt(self.filePointer, posn[np.newaxis,n,:] )
        #
        nbElements = 0
        maxDimensionalityOfelements = 0
        for ntype,elems in meshObject.elements.iteritems():
            maxDimensionalityOfelements = max(EN.dimension[ntype],maxDimensionalityOfelements)
            if EN.dimension[ntype] == maxDimensionalityOfelements:
                nbElements += elems.GetNumberOfElements()


        self.filePointer.write("  **element\n")
        self.filePointer.write("{}\n".format(nbElements))


        cpt =0;
        for ntype, data in meshObject.elements.iteritems():
            elemtype = GeofName[ntype]
            if EN.dimension[ntype] != maxDimensionalityOfelements and lowerDimElementsAsSets:
                continue
            #npe = data.GetNumberOfNodesPerElement()
            for i in xrange(data.GetNumberOfElements() ):
                if useOriginalId:
                    self.filePointer.write("{} {} ".format(data.originalIds[i],elemtype) )
                    self.filePointer.write(" ".join([str(int(meshObject.originalIDNodes[x])) for x in data.connectivity[i,:].flatten()]))
                else:
                    self.filePointer.write("{} {} ".format(cpt+1,elemtype) )
                    self.filePointer.write(" ".join([str(x+1) for x in data.connectivity[i,:].flatten()]))
                cpt += 1;
                self.filePointer.write("\n")

        self.filePointer.write(" ***group \n")

        for tag in meshObject.nodesTags:
            tag.tighten()
            self.filePointer.write("  **nset {} \n".format(tag.name))
            data = np.zeros((meshObject.GetNumberOfNodes(),1),dtype=np.int)
            if useOriginalId:
                self.filePointer.write(" ".join([str(meshObject.originalIDNodes[x]) for x in tag.id]))
            else:
                self.filePointer.write(" ".join([str(x+1) for x in tag.id]))
            self.filePointer.write("\n")

        meshObject.PrepareForOutput();

        if lowerDimElementsAsSets :
            celtags = meshObject.GetNamesOfCellTags()
            for tagname in celtags:

                idInTag = 0
                flag = False
                for ntype,elems in meshObject.elements.iteritems():
                    if EN.dimension[ntype] == maxDimensionalityOfelements:
                        if elems.tags.has_key(tagname):
                            flag =  True
                            idInTag += elems.tags[tagname].cpt
                self.PrintVerbose("Tag " + str(tagname) + " has "+ str(idInTag) + " elements")
                # no output if no elements in tag
                if flag == False : continue
                #empty tags
                if idInTag == 0 :  continue
                self.filePointer.write("  **elset {} \n".format(tagname))
                cpt =0

                for ntype,elems in meshObject.elements.iteritems():
                    if EN.dimension[ntype] != maxDimensionalityOfelements:
                        continue
                    if elems.tags.has_key(tagname):
                        tag = elems.tags[tagname]
                        tag.tighten()
                        if tag.cpt :
                            self.filePointer.write(" ".join([str(x+1+cpt) for x in tag.id]))
                    cpt += elems.GetNumberOfElements()

                self.filePointer.write("\n")
        else:
            celtags = meshObject.GetNamesOfCellTags()
            for tagname in celtags:
                self.filePointer.write("  **elset {} \n".format(tagname))
                data = meshObject.GetElementsInTag(tagname,useOriginalId=useOriginalId)
                if useOriginalId :
                    self.filePointer.write(" ".join([str(x) for x in data]))
                else:
                    self.filePointer.write(" ".join([str(x+1) for x in data]))
                self.filePointer.write("\n")

        # Dotsets, lisets, facets
        if lowerDimElementsAsSets:
            for dimToTreat in xrange(maxDimensionalityOfelements):

                celtags = meshObject.GetNamesOfCellTags()
                for tagname in celtags:

                    idInTag = 0
                    flag = False
                    for ntype,elems in meshObject.elements.iteritems():
                        if EN.dimension[ntype] == dimToTreat:
                            if elems.tags.has_key(tagname):
                                flag =  True
                                idInTag += elems.tags[tagname].cpt
                    self.PrintVerbose("Set  " + str(tagname) + " has "+ str(idInTag) + " elements")
                    # no output if no elements in tag
                    if flag == False : continue
                    #empty tags
                    if idInTag == 0 :  continue
                    if dimToTreat == 0:
                        self.filePointer.write("  **doset {} \n".format(tagname))
                    elif dimToTreat == 1:
                        self.filePointer.write("  **liset {} \n".format(tagname))
                    elif dimToTreat == 2:
                        self.filePointer.write("  **faset {} \n".format(tagname))


                    for ntype,elems in meshObject.elements.iteritems():
                        if EN.dimension[ntype] != dimToTreat:
                            continue
                        if elems.tags.has_key(tagname):
                            tag = elems.tags[tagname]
                            tag.tighten()
                            name = GeofSetName[ntype];


                            for e in xrange(tag.cpt):
                                self.filePointer.write(" {} ".format(name))
                                self.filePointer.write(" ".join([str(x+1) for x in elems.connectivity[tag.id[e],:] ]))
                                self.filePointer.write(" \n")


                    self.filePointer.write("\n")


        self.filePointer.write("****return \n")

def CheckIntegrity():
    import OTTools.FE.UnstructuredMesh as UM

    from OTTools.Helpers.Tests import TestTempDir

    tempdir = TestTempDir.GetTempPath()

    mymesh = UM.UnstructuredMesh()
    mymesh.nodes = np.array([[0.00000000001,0,0],[1,0,0],[0,1,0],[1,1,0]],dtype=np.float)
    mymesh.originalIDNodes = np.array([1, 3, 4, 5],dtype=np.int)

    mymesh.nodesTags.CreateTag("coucou").AddToTag(0)

    tris = mymesh.GetElementsOfType(EN.Triangle_3)
    tris.AddNewElement([0,1,2],0)
    tris.AddNewElement([2,1,3],3)
    tris.originalIds = np.array([3, 5],dtype=np.int)

    tris = mymesh.GetElementsOfType(EN.Bar_2)
    tris.AddNewElement([0,1],0)
    tris.AddNewElement([1,3],1)


    mymesh.AddElementToTagUsingOriginalId(3,"Tag1")
    mymesh.AddElementToTagUsingOriginalId(5,"Tag3")

    OW = GeofWriter()
    print(OW)
    OW.Open(tempdir+"Test_GeoWriter.geof")
    OW.Write(mymesh, useOriginalId=True)
    OW.Close()

    WriteMeshToGeof(tempdir+"Test_GeoWriter_II.geof", mymesh,useOriginalId=False )
    WriteMeshToGeof(tempdir+"Test_GeoWriter_III.geof", mymesh,useOriginalId=False,lowerDimElementsAsSets=True )
    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
