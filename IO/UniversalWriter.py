# -*- coding: utf-8 -*-

__author__ = "Felipe Bordeu"

externalWriters = {}

def WriterFactory(nameOrFilename,ops={"default":"xmf"}):

    import os.path
    dirname = os.path.dirname(nameOrFilename)
    basename,extention = os.path.splitext(os.path.basename(nameOrFilename))

    res = None
    if extention == ".geof" or nameOrFilename == "geof":
        from BasicTools.IO.GeofWriter import GeofWriter
        res = GeofWriter()
        res.SetWriteLowerDimElementsAsSets(True)
    elif extention ==  ".gmsh" or nameOrFilename == "gmsh":
        from  BasicTools.IO.GmshWriter import GmshWriter
        res = GmshWriter()
    elif extention ==  ".mesh" or nameOrFilename == "mesh":
        from BasicTools.IO.MeshWriter import MeshWriter
        res = MeshWriter()
    elif extention ==  ".xmf" or extention ==  ".xdmf" or nameOrFilename == "xmf":
        from BasicTools.IO.XdmfWriter import  XdmfWriter
        res = XdmfWriter()
    else:
        if extention in externalWriters:
            res = externalWriters[extention]()
        else:
            res = WriterFactory(ops["default"])

    res.SetFileName(nameOrFilename)
    return res


def WriteMesh(filename,out,binary=False):# pragma: no cover

    writer = WriterFactory(filename)
    writer.SetBinary(binary)
    writer.Open()
    writer.write(out)
    writer.Close()

## to use this function add this lines to the
## programmmable filter in paraview, and change the output type to UnstructuredGrid
##
#
#filename = "here you put your filename"
#from BasicTools.IO.UniversalWriter import PopulateMeshFromVtkAndWriteMesh
#PopulateMeshFromVtkAndWriteMesh(filename,self.GetInput())

def PopulateMeshFromVtkAndWriteMesh(filename, vtkobject):# pragma: no cover

    from BasicTools.FE.UnstructuredMeshTools import VtkToMesh
    mesh = VtkToMesh(vtkobject)
    WriteMesh(filename,mesh)

def CheckIntegrity():
    print(WriterFactory("toto.geof"))
    print(WriterFactory("toto.gmsh"))
    print(WriterFactory("toto.mesh"))
    print(WriterFactory("toto.xdmf"))
    print(WriterFactory("toto"))
    #normally this class must have the same API as
    #
    from BasicTools.IO.WriterBase import WriterBase as WriterBase

    class MyCustomWriter(WriterBase):
        pass

    externalWriters[".myExtention"] = MyCustomWriter
    print(WriterFactory("toto.myExtention"))

    return "ok"


if __name__ == '__main__':
    print((CheckIntegrity()))# pragma: no cover