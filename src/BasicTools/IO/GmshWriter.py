# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

""" Gmsh  File Writer (gmesh mesh files)

"""
import numpy as np
import BasicTools.Containers.ElementNames as EN

from BasicTools.Containers.MeshBase import Tag as Tag
from BasicTools.IO.WriterBase import WriterBase as WriterBase
from BasicTools.NumpyDefs import PBasicFloatType, PBasicIndexType
from BasicTools.IO.GmshTools import gmshName,PermutationGmshToBasicTools
PermutationBasicToolsToGmsh = {key:np.argsort(value) for key, value in PermutationGmshToBasicTools.items() }

def WriteMeshToGmsh(filename,mesh, useOriginalId=False):
    OW = GmshWriter()
    OW.Open(filename)
    OW.Write(mesh,useOriginalId = useOriginalId)
    OW.Close()

class GmshWriter(WriterBase):
    def __init__(self):
        super(GmshWriter,self).__init__()
        self.SetBinary(False)
        self.canHandleBinaryChange = False

    def __str__(self):
        res  = 'GmshWriter : \n'
        res += '   FileName : '+str(self.fileName)+'\n'
        return res

    def SetFileName(self,fileName):
        self.fileName = fileName;

    def Write(self,meshObject,useOriginalId=False,PointFieldsNames=None,PointFields=None,CellFieldsNames=None,CellFields=None):
        self.filePointer.write("$MeshFormat\n");
        self.filePointer.write("2.2 0 8\n");
        self.filePointer.write("$EndMeshFormat\n");
        self.filePointer.write("$Nodes\n");
        numberofpoints = meshObject.GetNumberOfNodes()
        self.filePointer.write("{} \n".format(numberofpoints) )
        #
        posn = meshObject.GetPosOfNodes()
        if useOriginalId:
           for n in range(numberofpoints):
               self.filePointer.write("{} ".format(int(meshObject.originalIDNodes[n])))
               #np.savetxt(self.filePointer, posn[np.newaxis,n,:] )
               posn[np.newaxis,n,:].tofile(self.filePointer, sep=" ")
               self.filePointer.write("\n")
        else:
           for n in range(numberofpoints):
               self.filePointer.write("{} ".format(n+1) )
               #np.savetxt(self.filePointer, posn[np.newaxis,n,:] )
               posn[np.newaxis,n,:].tofile(self.filePointer, sep=" ")
               self.filePointer.write("\n")
        self.filePointer.write("$EndNodes\n");
        self.filePointer.write("$Elements\n");
        self.filePointer.write("{}\n".format(meshObject.GetNumberOfElements()))

        cpt = 1

        #meshObject.PrepareForOutput();

        celltags = meshObject.GetNamesOfElemTags()
        elements = meshObject.elements
        tagcounter = 2
        #for tagname in celtags:
        for elementContainer in elements:
            elemtype = gmshName[elementContainer]
            #try:
            data = meshObject.elements[elementContainer]
            if useOriginalId:
                 for i in range(data.GetNumberOfElements() ):
                    self.filePointer.write("{}{} {} {} {} ".format(data.originalIds[i],elemtype,2,1,1) )
                    self.filePointer.write(" ".join([str(meshObject.originalIDNodes[x]) for x in data.connectivity[i,:].ravel()]))
                    cpt += 1
                    self.filePointer.write("\n")
            else:
                #for connectivity in meshObject.elements[elementContainer].connectivity[meshObject.elements[elementContainer].tags[tagname].id-1]:
                for i in range(data.GetNumberOfElements() ):
                  connectivity = data.connectivity[i,:]
                  if elementContainer in PermutationBasicToolsToGmsh:
                      connectivity = [connectivity[x] for x in PermutationBasicToolsToGmsh[elementContainer]]
                  self.filePointer.write("{} {} {} {} {} ".format(cpt,elemtype,2,tagcounter,tagcounter) )
                  self.filePointer.write(" ".join([str(x+1) for x in connectivity]))
                  cpt += 1
                  self.filePointer.write("\n")
            #except KeyError:
            #  continue
          #tagcounter += 1

        self.filePointer.write("$EndElements\n")

from BasicTools.IO.IOFactory import RegisterWriterClass
RegisterWriterClass(".msh",GmshWriter)

def CheckIntegrity():
    import BasicTools.Containers.UnstructuredMesh as UM

    from BasicTools.Helpers.Tests import TestTempDir

    tempdir = TestTempDir.GetTempPath()

    mymesh = UM.UnstructuredMesh()
    mymesh.nodes = np.array([[0.00000000001,0,0],[1,0,0],[0,1,0],[1,1,0]],dtype=PBasicFloatType)
    mymesh.originalIDNodes = np.array([1, 3, 4, 5],dtype=PBasicIndexType)


    mymesh.nodesTags.CreateTag("coucou").AddToTag(0)

    tris = mymesh.GetElementsOfType('tri3')
    tris.AddNewElement([0,1,2],0)
    tris.AddNewElement([2,1,3],3)
    tris.originalIds = np.array([3, 5],dtype=PBasicIndexType)

    mymesh.AddElementToTagUsingOriginalId(3,"Tag1")
    mymesh.AddElementToTagUsingOriginalId(5,"Tag3")

    OW = GmshWriter()
    print(OW)
    OW.Open(tempdir+"Test_GmshWriter.geof")
    OW.Write(mymesh, useOriginalId=False)
    OW.Close()


    print(mymesh)
    WriteMeshToGmsh(tempdir+"Test_GmshWriter_II.geof", mymesh,useOriginalId=True  )
    #print(tempdir)
    #print(open(tempdir+"Test_GmshWriter_II.geof").read())

    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
