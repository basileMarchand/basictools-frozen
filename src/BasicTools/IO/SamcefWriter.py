# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

import numpy as np

from BasicTools.IO.WriterBase import WriterBase as WriterBase
import BasicTools.Containers.ElementNames as EN

def WriteMeshToDat(filename,mesh, normals= None):
    OW = DatWriter()
    OW.Open(filename)
    OW.Write(mesh, normals=normals)
    OW.Close()

BasicToolToSamcef = dict()
BasicToolToSamcef[EN.Bar_2] = (np.array([1,1]), np.arange(2),2)
BasicToolToSamcef[EN.Triangle_3] = (np.array([1,1,1]), np.arange(3),3)
BasicToolToSamcef[EN.Tetrahedron_4] = (np.array([1,1,1,1]), np.arange(4),4)

class DatWriter(WriterBase):
    def __init__(self):
        super(DatWriter,self).__init__()
    def __str__(self):
        res  = 'DatWriter : \n'
        res += '   FileName : '+str(self.fileName)+'\n'
        return res

    def Write(self,meshObject,normals=None,Name= None,PointFieldsNames=None,PointFields=None,CellFieldsNames=None,CellFields=None,GridFieldsNames=None,GridFields=None):
        #Nodes
        numberofpoints = meshObject.GetNumberOfNodes()
        posn = meshObject.GetPosOfNodes()
        self.writeText(".NOE\n")
        for n in range(numberofpoints):
            self.filePointer.write("{} ".format(n+1))
            posn[np.newaxis,n,:].tofile(self.filePointer, sep=" ")
            self.writeText("\n")

        groupcpt = 1
        gropuNames = dict()
        #Nodes Tags
        for tag in meshObject.nodesTags:
            self.writeText('.SEL GROUP {} NOM "{}" NOEUDS \n'.format(groupcpt,tag.name) )
            gropuNames[tag.name] =  groupcpt
            groupcpt += 1
            self.writeText('I ')
            #tag meshObject.nodesTags[tagname]
            (tag.GetIds()+1).tofile(self.filePointer, sep=" ")
            self.writeText("\n")

        #Elements
        self.writeText(".MAI\n")
        cpt =0
        #for tagname in celtags:
        for name,data in meshObject.elements.items():
            sign,permutation,splitpoint = BasicToolToSamcef[name]

            lconn = (1+data.connectivity)*sign
            lconn = lconn[permutation]

            fp = lconn[:,:splitpoint]
            sp = lconn[:,splitpoint:]
            for n in range(data.GetNumberOfElements()):
                self.writeText("I {} N ".format(cpt))
                lconn = (1+data.connectivity[n,:])*sign

                fp[n,:].tofile(self.filePointer, sep=" ")
                if sp.shape[1] > 0:
                    self.writeText(" 0 ")
                    sp[n,:].tofile(self.filePointer, sep=" ")
                cpt += 1
                self.writeText("\n")

        celltags = meshObject.GetNamesOfElemTags()
        #Nodes Tags
        for tagname in celltags:
            self.writeText('.SEL GROUP {} NOM "{}" MAILLE \n'.format(groupcpt,tagname) )
            gropuNames[tag.name] =  groupcpt
            groupcpt += 1

            self.writeText('I ')
            ids = meshObject.GetElementsInTag(tagname)+1
            ids.tofile(self.filePointer, sep=" ")
            self.writeText("\n")

from BasicTools.IO.IOFactory import RegisterWriterClass
RegisterWriterClass(".datt",DatWriter)

def CheckIntegrity():
    import BasicTools.Containers.UnstructuredMesh as UM

    from BasicTools.Helpers.Tests import TestTempDir

    tempdir = TestTempDir.GetTempPath()

    mymesh = UM.UnstructuredMesh()
    mymesh.nodes = np.array([[0.00000000001,0,0],[1,0,0],[0,1,0],[1,1,0]],dtype=np.float)
    mymesh.originalIDNodes = np.array([1, 3, 4, 5],dtype=np.int)

    mymesh.nodesTags.CreateTag("FirstNode").AddToTag(0)

    tris = mymesh.GetElementsOfType(EN.Triangle_3)
    tris.AddNewElement([0,1,2],0)
    tris.AddNewElement([2,1,3],3)
    tris.originalIds = np.array([3, 5],dtype=np.int)

    mymesh.AddElementToTagUsingOriginalId(3,"Tag1")
    mymesh.AddElementToTagUsingOriginalId(5,"Tag3")

    OW = DatWriter()
    OW.Open(tempdir+"Test_GmshWriter.datt")
    OW.Write(mymesh)
    OW.Close()

    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
