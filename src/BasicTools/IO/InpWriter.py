# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

"""Inp file writer (Abaqus mesh files)
"""

import numpy as np

import BasicTools.Containers.ElementNames as EN
from BasicTools.IO.WriterBase import WriterBase as WriterBase
from BasicTools.IO.AbaqusTools import InpNameToBasicTools, permutation, BasicToolsToInpName

def WriteMeshToINP(filename,mesh, useOriginalId=False):
    """Export Mesh to disk in the inp (Abaqus mesh) format file.
        A file is created using the path and name of filename

    Parameters
    ----------
    filename : str
        name with path to the file to be created (relative or absolute)
    mesh : UnstructuredMesh
        the mesh to be exported
    useOriginalId : bool, optional
        If True, Original Id for the number of nodes and elements are used
        (the user is responsible of the consistency of this data), by default False
    """
    OW = InpWriter()
    OW.Open(filename)
    OW.Write(mesh,useOriginalId = useOriginalId)
    OW.Close()

class InpWriter(WriterBase):
    """Class to writes a inp file on disk
    """
    def __init__(self):
        super(InpWriter,self).__init__()
        self.canHandleBinaryChange = False
        self.canHandleTemporal = False
        self.canHandleMultidomain = False

    def __str__(self):
        res  = 'InpWriter : \n'
        res += '   FileName : '+ str(self.fileName) +'\n'
        return res

    def Write(self, meshObject, useOriginalId=False, PointFieldsNames=None, PointFields=None, CellFieldsNames=None, CellFields=None):
        """Function to writes a CGNS File on disk

        Parameters
        ----------
        mesh : UnstructuredMesh
            the mesh to be exported
        useOriginalId : bool, optional
            If True, Original Id for the number of nodes and elements are used
            (the user is responsible of the consistency of this data), by default False
        PointFieldsNames : None
            Not Used, by default None
        PointFields : None
            Not Used, by default None
        CellFieldsNames : None
            Not Used, by default None
        CellFields : None
            Not Used, by default None
        """
        if PointFieldsNames is not None or \
            PointFields      is not None or \
            CellFieldsNames  is not None or \
            CellFields       is not None:
                print("warning InpWriter only can write the mesh, fields are ignored")

        meshObject.PrepareForOutput()

        self.filePointer.write("** Written by BasicTools package\n")
        self.filePointer.write("*NODE\n");

        numberofpoints = meshObject.GetNumberOfNodes()
        #self.filePointer.write("{} {} \n".format(numberofpoints,meshObject.GetDimensionality()) )
        #
        posn = meshObject.GetPosOfNodes()
        if useOriginalId:
            for n in range(numberofpoints):
                self.filePointer.write("{}, ".format(int(meshObject.originalIDNodes[n])))
                #np.savetxt(self.filePointer, posn[n,:] )
                posn[n,:].tofile(self.filePointer, sep=", ")
                self.filePointer.write("\n")
        else:
            for n in range(numberofpoints):
                self.filePointer.write("{}, ".format(n+1) )
                #np.savetxt(self.filePointer, posn[np.newaxis,n,:] )
                posn[np.newaxis,n,:].tofile(self.filePointer, sep=", ")
                self.filePointer.write("\n")
        #


        cpt =0
        FENames = None
        if "FE Names" in meshObject.elemFields.keys():
            FENames = meshObject.elemFields["FE Names"]

        for elemtype, data in meshObject.elements.items():

            if FENames is None:
                self.filePointer.write(f"*ELEMENT, TYPE={BasicToolsToInpName[elemtype]}\n")
            else:
                if np.any(FENames[cpt] != FENames[cpt:data.GetNumberOfElements()] ):
                    raise(Exception("Error, heterogeneous FE Names not supported yet sorry!!"))
                self.filePointer.write(f"*ELEMENT, TYPE={FENames[cpt]}\n")

            for i in range(data.GetNumberOfElements() ):
                conn = data.connectivity[i,:].ravel()
                if useOriginalId:
                    self.filePointer.write("{}, ".format(data.originalIds[i]) )
                    self.filePointer.write(", ".join([str(meshObject.originalIDNodes[x]  ) for x in conn]))
                else:
                    self.filePointer.write("{}, ".format(cpt+1) )
                    self.filePointer.write(", ".join(map(str,conn+1)))
                cpt += 1
                self.filePointer.write("\n")

        for tag in meshObject.nodesTags:
            if len(tag) == 0:
                continue
            self.filePointer.write("*NSET, NSET={} \n".format(tag.name))
            if useOriginalId:
                self.filePointer.write(", ".join([str(int(meshObject.originalIDNodes[x])) for x in tag.GetIds()]))
            else:
                self.filePointer.write(", ".join([str(x+1) for x in tag.GetIds()]))
            self.filePointer.write("\n")

        elemtags = meshObject.GetNamesOfElemTags()
        for tagname in elemtags:
            self.filePointer.write("*ELSET, ELSET={} \n".format(tagname))
            data = meshObject.GetElementsInTag(tagname,useOriginalId=useOriginalId)
            if useOriginalId :
                self.filePointer.write(" ".join([str(x) for x in data]))
                self.filePointer.write(" ")
            else:
                self.filePointer.write(" ".join([str(x+1) for x in data]))
                self.filePointer.write(" ")
            self.filePointer.write("\n")


from BasicTools.IO.IOFactory import RegisterWriterClass
RegisterWriterClass(".inp",InpWriter)

def CheckIntegrity():
    from BasicTools.Helpers.Tests import TestTempDir
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateCube

    tempdir = TestTempDir.GetTempPath()
    print(tempdir)

    mesh = CreateCube()
    mesh.GenerateManufacturedOriginalIDs()


    head = "** this is the head of the file \n"
    tail = "** and this is the tail of this file\n"
    mesh.nodes +=1
    ids = [0]
    mesh.nodes[ids,0] *= 2

    OW = InpWriter()
    OW.Open(tempdir+"Test_InpWriter.inp")
    OW.writeText(head)
    OW.Write(mesh, useOriginalId=True)
    OW.writeText(tail)
    OW.Close()
    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
