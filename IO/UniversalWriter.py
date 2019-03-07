# -*- coding: utf-8 -*-

__author__ = "Felipe Bordeu"

#externalWriters = {}
#
#def WriterFactory(nameOrFilename,ops={"default":"xmf"}):
#
#    import os.path
#    dirname = os.path.dirname(nameOrFilename)
#    basename,extention = os.path.splitext(os.path.basename(nameOrFilename))
#
#    res = None
#    if extention == ".geof" or nameOrFilename == "geof":
#        from BasicTools.IO.GeofWriter import GeofWriter
#        res = GeofWriter()
#        res.SetWriteLowerDimElementsAsSets(True)
#    elif extention ==  ".gmsh" or nameOrFilename == "gmsh":
#        from  BasicTools.IO.GmshWriter import GmshWriter
#        res = GmshWriter()
#    elif extention ==  ".mesh" or nameOrFilename == "mesh":
#        from BasicTools.IO.MeshWriter import MeshWriter
#        res = MeshWriter()
#    elif extention ==  ".xmf" or extention ==  ".xdmf" or nameOrFilename == "xmf":
#        from BasicTools.IO.XdmfWriter import XdmfWriter
#        res = XdmfWriter()
#
#    elif extention ==  ".ut" :
#
#        import BasicTools.IO.UtWriter as UW
#        res = UW.UtWriter()
#        import numpy as np
#        res.AttachSequence(np.array([[0,0,0,0,0]]))
#    else:
#        if extention in externalWriters:
#            res = externalWriters[extention]()
#        else:
#            raise("Unable to find a suitable writer for the file  :"+str() )
#            #res = WriterFactory(ops["default"])
#
#    res.SetFileName(nameOrFilename)
#    return res


def WriteMesh(filename,outmesh,binary=False):# pragma: no cover

    from BasicTools.IO.IOFactory import CreateWriter
    writer = CreateWriter("."+filename.split(".")[-1])
    writer.SetFileName(filename)
    writer.SetBinary(binary)
    writer.Open()

    PointFields = None
    PointFieldsNames = None
    if hasattr(outmesh,"nodeFields"):
        PointFieldsNames = list(outmesh.nodeFields.keys())
        PointFields = list(outmesh.nodeFields.values())

    CellFields = None
    CellFieldsNames = None
    if hasattr(outmesh,"elemFields"):
        CellFieldsNames = list(outmesh.elemFields.keys())
        CellFields = list(outmesh.elemFields.values())

    writer.Write(outmesh,PointFieldsNames=PointFieldsNames,PointFields=PointFields,CellFieldsNames=CellFieldsNames,CellFields=CellFields )
    writer.Close()

## to use this function add this lines to the
## programmmable filter in paraview, and change the output type to UnstructuredGrid
##
#
#filename = "here you put your filename"
#from BasicTools.IO.UniversalWriter import PopulateMeshFromVtkAndWriteMesh
#PopulateMeshFromVtkAndWriteMesh(filename,self.GetInput())

def PopulateMeshFromVtkAndWriteMesh(filename, vtkobject):# pragma: no cover

    from BasicTools.Containers.vtkBridge import VtkToMesh
    mesh = VtkToMesh(vtkobject)

    WriteMesh(filename,mesh)

def CheckIntegrity():
    from BasicTools.IO.IOFactory import CreateWriter, InitAllWriters, RegisterWriterClass
    InitAllWriters()

    print(CreateWriter(".geof"))
    print(CreateWriter(".msh"))
    print(CreateWriter(".mesh"))
    print(CreateWriter(".xdmf"))
    try:
        print(CreateWriter("toto"))
    except:
        pass
    else:
        raise (Exception("this must fail " ))

    #normally this class must have the same API as
    #
    from BasicTools.IO.WriterBase import WriterBase as WriterBase

    class MyCustomWriter(WriterBase):
        pass

    RegisterWriterClass(".myExtention",MyCustomWriter)

    print(CreateWriter(".myExtention"))

    return "ok"


if __name__ == '__main__':
    print((CheckIntegrity()))# pragma: no cover
