# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

import numpy as np

from BasicTools.IO.WriterBase import WriterBase as WriterBase
import BasicTools.Containers.ElementNames as EN
from BasicTools.NumpyDefs import PBasicIndexType
from BasicTools.NumpyDefs import PBasicFloatType

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

    def Write(self,meshObject,PointFieldsNames=None,PointFields=None,CellFieldsNames=None,CellFields=None,GridFieldsNames=None,GridFields=None):
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
        if len(meshObject.nodesTags):
            self.writeText('.SEL\n' )
            for tag in meshObject.nodesTags:
                self.writeText('GROUP {} NOM "{}" NOEUDS \n'.format(groupcpt,tag.name) )
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
            lconn = lconn[:,permutation]

            fp = lconn[:,:splitpoint]
            sp = lconn[:,splitpoint:]
            for n in range(data.GetNumberOfElements()):
                self.writeText("I {} N ".format(cpt))
                fp[n,:].tofile(self.filePointer, sep=" ")
                if sp.shape[1] > 0:
                    self.writeText(" 0 ")
                    sp[n,:].tofile(self.filePointer, sep=" ")
                cpt += 1
                self.writeText("\n")

        celltags = meshObject.GetNamesOfElemTags()
        #Nodes Tags
        if len(celltags):
            self.writeText('.SEL\n' )
            for tagname in celltags:
                self.writeText('GROUP {} NOM "{}" MAILLE \n'.format(groupcpt,tagname) )
                gropuNames[tagname] =  groupcpt
                groupcpt += 1

                self.writeText('I ')
                ids = meshObject.GetElementsInTag(tagname)+1
                ids.tofile(self.filePointer, sep=" ")
                self.writeText("\n")

        if "FrameX_0" in CellFieldsNames:
            FrameX_0 = CellFields[CellFieldsNames.index("FrameX_0")]
            FrameX_1 = CellFields[CellFieldsNames.index("FrameX_1")]
            FrameX_2 = CellFields[CellFieldsNames.index("FrameX_2")]
            FrameY_0 = CellFields[CellFieldsNames.index("FrameY_0")]
            FrameY_1 = CellFields[CellFieldsNames.index("FrameY_1")]
            FrameY_2 = CellFields[CellFieldsNames.index("FrameY_2")]

            nbelems = meshObject.GetNumberOfElements()

            self.writeText('.FRA\n')
            idx = []
            for n in range(nbelems):
                if np.sum(np.abs([FrameX_0[n],FrameX_1[n],FrameX_2[n] ])) == 0:
                    continue
                idx.append(n+1)
                self.writeText(f'I {n+1} V1 {FrameX_0[n]} {FrameX_1[n]} {FrameX_2[n]} V2 {FrameY_0[n]} {FrameY_1[n]} {FrameY_2[n]} \n')

            self.writeText('.AEL\n')
            for n in idx:
                self.writeText(f'I {n} FRAME {n} \n')
        self.writeText('RETURN\n')



from BasicTools.IO.IOFactory import RegisterWriterClass
RegisterWriterClass(".datt",DatWriter)

def CheckIntegrity():
    import BasicTools.Containers.UnstructuredMesh as UM

    from BasicTools.Helpers.Tests import TestTempDir

    tempdir = TestTempDir.GetTempPath()

    mymesh = UM.UnstructuredMesh()
    mymesh.nodes = np.array([[0.00000000001,0,0],[1,0,0],[0,1,0],[1,1,0]],dtype=PBasicFloatType)
    mymesh.originalIDNodes = np.array([1, 3, 4, 5],dtype=PBasicIndexType)

    mymesh.nodesTags.CreateTag("FirstNode").AddToTag(0)

    tris = mymesh.GetElementsOfType(EN.Triangle_3)
    tris.AddNewElement([0,1,2],0)
    tris.AddNewElement([2,1,3],3)
    tris.originalIds = np.array([3, 5],dtype=PBasicIndexType)

    mymesh.AddElementToTagUsingOriginalId(3,"Tag1")
    mymesh.AddElementToTagUsingOriginalId(5,"Tag3")

    nbel = mymesh.GetNumberOfElements()
    CellFields= [np.ones(nbel),np.zeros(nbel),np.zeros(nbel),np.zeros(nbel),np.ones(nbel),np.zeros(nbel)]
    CellFieldsNames= ["FrameX_0","FrameX_1","FrameX_2","FrameY_0","FrameY_1","FrameY_2"]

    OW = DatWriter()
    OW.Open(tempdir+"Test_SamcefWriter.datt")
    OW.Write(mymesh,CellFieldsNames=CellFieldsNames,CellFields=CellFields)
    OW.Close()

    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
