#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
# -*- coding: utf-8 -*-
import numpy as np

from BasicTools.IO.WriterBase import WriterBase as WriterBase
import BasicTools.Containers.ElementNames as EN


def WriteMeshToStl(filename,mesh, normals= None):
    OW = StlWriter()
    OW.Open(filename)
    OW.Write(mesh, normals=normals)
    OW.Close()

class StlWriter(WriterBase):
    def __init__(self):
        super(StlWriter,self).__init__()
        self.extractSkin = True

    def __str__(self):
        res  = 'StlWriter : \n'
        res += '   FileName : '+str(self.fileName)+'\n'
        return res

    def SetFileName(self,fileName):
        self.fileName = fileName;

    def Close(self):
        if self.isOpen():
            super(StlWriter,self).Close()

    def Write(self,meshObject,normals=None,Name= None,PointFieldsNames=None,PointFields=None,CellFieldsNames=None,CellFields=None,GridFieldsNames=None,GridFields=None):

        if self.extractSkin :
            from BasicTools.Containers.UnstructuredMeshTools import ComputeSkin, CleanLonelyNodes
            skin = ComputeSkin(meshObject)
            #print(skin)
            CleanLonelyNodes(skin)
            nodes = skin.nodes
            elements = skin.GetElementsOfType(EN.Triangle_3)
        else:
            nodes = meshObject.nodes
            elements = meshObject.GetElementsOfType(EN.Triangle_3)
        #self.PrintDebug(Name)
        if Name is None:
            Name == "Solid"

        self.writeText("solid {}\n".format(Name))
        for i in range(elements.GetNumberOfElements()):
            if normals is not None:
                self.writeText(" facet normal {}\n".format(" ".join(map(str,normals[i,:])) ))
            else:
                p0 = nodes[elements.connectivity[i,0],:]
                p1 = nodes[elements.connectivity[i,1],:]
                p2 = nodes[elements.connectivity[i,2],:]
                normal = np.cross(p1-p0,p2-p0)
                normal = normal/np.linalg.norm(normal)
                self.writeText(" facet normal {}\n".format(" ".join(map(str,normal)) ))
            self.writeText("  outer loop\n")
            for p in range(3):
                self.writeText("   vertex {}\n".format(" ".join(map(str,nodes[elements.connectivity[i,p],:] ))))
            self.writeText("  endloop\n")
            self.writeText(" endfacet\n")
        self.writeText("endsolid\n")

from BasicTools.IO.IOFactory import RegisterWriterClass
RegisterWriterClass(".stl",StlWriter)


def CheckIntegrity():
    data = u"""   solid cube_corner
          facet normal 0.0 -1.0 0.0
            outer loop
              vertex 0.0 0.0 0.0
              vertex 1.0 0.0 0.0
              vertex 0.0 0.0 1.0
            endloop
          endfacet
          facet normal 0.0 0.0 -1.0
            outer loop
              vertex 0.0 0.0 0.0
              vertex 0.0 1.0 0.0
              vertex 1.0 0.0 0.0
            endloop
          endfacet
          facet normal -1.0 0.0 0.0
            outer loop
              vertex 0.0 0.0 0.0
              vertex 0.0 0.0 1.0
              vertex 0.0 1.0 0.0
            endloop
          endfacet
          facet normal 0.577 0.577 0.577
            outer loop
              vertex 1.0 0.0 0.0
              vertex 0.0 1.0 0.0
              vertex 0.0 0.0 1.0
            endloop
          endfacet
        endsolid"""
    from BasicTools.IO.StlReader import ReadStl as ReadStl
    from BasicTools.Helpers.Tests import TestTempDir

    res = ReadStl(string=data)
    print(res)
    tempdir = TestTempDir.GetTempPath()
    WriteMeshToStl(tempdir+"Test_StlWriter.stl",res)
    WriteMeshToStl(tempdir+"Test_StlWriter_with_normals.stl",res, normals =res.elemFields["normals"])

    ow = StlWriter()
    print(ow)

    return("ok")

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
