# -*- coding: utf-8 -*-
from BasicTools.IO.WriterBase import WriterBase as WriterBase
import BasicTools.FE.ElementNames as EN
import numpy as np

def WriteMeshToStl(filename,mesh, normals= None):
    OW = StlWriter()
    OW.Open(filename)
    OW.Write(mesh, normals=normals)
    OW.Close()

class StlWriter(WriterBase):
    def __init__(self):
        super(StlWriter,self).__init__()

    def __str__(self):
        res  = 'StlWriter : \n'
        res += '   FileName : '+str(self.fileName)+'\n'
        return res

    def SetFileName(self,fileName):
        self.fileName = fileName;

    def Write(self,meshObject,normals=None,Name= None):

        elements = meshObject.GetElementsOfType(EN.Triangle_3)

        self.filePointer.write("solid {}\n".format(Name))
        for i in xrange(elements.GetNumberOfElements()):
            if normals is not None:
                self.filePointer.write(" facet normal {}\n".format(" ".join(map(str,normals[i,:])) ))
            else:
                p0 = meshObject.nodes[elements.connectivity[i,0],:]
                p1 = meshObject.nodes[elements.connectivity[i,1],:]
                p2 = meshObject.nodes[elements.connectivity[i,2],:]
                normal = np.cross(p1-p0,p2-p0)
                normal = normal/np.linalg.norm(normal)
                self.filePointer.write(" facet normal {}\n".format(" ".join(map(str,normal)) ))
            self.filePointer.write("  outer loop\n")
            for p in xrange(3):
                self.filePointer.write("   vertex {}\n".format(" ".join(map(str,meshObject.nodes[elements.connectivity[i,p],:] ))))
            self.filePointer.write("  endloop\n")
            self.filePointer.write(" endfacet\n")
        self.filePointer.write("endsolid\n")


def CheckIntegrity():
    data = """   solid cube_corner
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

    res,normals = ReadStl(string=data)
    tempdir = TestTempDir.GetTempPath()
    WriteMeshToStl(tempdir+"Test_StlWriter.stl",res)
    WriteMeshToStl(tempdir+"Test_StlWriter_with_normals.stl",res, normals=normals)

    ow = StlWriter()
    print(ow)

    return("ok")

if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
