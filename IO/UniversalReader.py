# -*- coding: utf-8 -*-

import OTTools.FE.ElementNames as ElementNames




def ReadMesh(filename):
    extention = filename.split(".")[-1].lower()
    if extention ==  "asc":
        import OTTools.IO.AscReader as AscReader
        return AscReader.ReadAsc(filename)
    elif extention ==  "geof":
        import OTTools.IO.GeofReader as GeofReader
        return GeofReader.ReadGeof(filename)
    elif extention ==  "msh":
        import OTTools.IO.GmshReader as GmshReader
        return GmshReader.ReadGmsh(filename)
    elif extention ==  "inp":
        import OTTools.IO.InpReader as ImpReader
        return ImpReader.ReadInp(filename)
    elif extention ==  "mesh":
        import OTTools.IO.MeshReader as ReadMesh
        return ReadMesh.ReadMesh(fileName=filename)
    elif extention ==  "meshb":
        import OTTools.IO.MeshReader as MeshReader
        return MeshReader.ReadMesh(fileName=filename)
    elif extention ==  "gcode":
        import OTTools.IO.GReader as GReader
        return GReader.ReadGCode(fileName=filename)
    elif extention ==  "solb":
        import OTTools.IO.MeshReader as MeshReader
        f = ".".join(filename.split(".")[0:-1]) + ".meshb"
        reader = MeshReader.MeshReader()
        reader.SetFileName(fileName=f)
        reader.Read()
        mesh = reader.output
        mesh.nodeFields = reader.ReadExtraFieldBinary(filename);
        return mesh
    else:
        raise Exception ("Unkown file extention : " + str(extention))


## to use this function add this lines to the
## programmmable source in paraview, and change the output type to UnstructuredGrid
##
#
#UniversalReader.ReadMeshAnPopulateVtkObject(filename,self.GetOutput())
#from OTTools.IO.UniversalReader import ReadMeshAndPopulateVtkObject as ReadMeshAndPopulateVtkObject
#filename = "here you put your filename"
#ReadMeshAndPopulateVtkObject(filename,self.GetOutput())

def ReadMeshAndPopulateVtkObject(filename, vtkobject= None):
    mesh = ReadMesh(filename)
    from OTTools.FE.UnstructuredTools import MeshToVtk
    return MeshToVtk(mesh, vtkobject)




def CheckIntegrity():
    from OTTools.FE.UnstructuredMeshTools import CreateMeshOfTriangles
    import numpy as np

    res = CreateMeshOfTriangles([[0,0,0],[1,0,0],[0,1,0],[0,0,1] ], [[0,1,2],[0,2,3]])
    res.elemFields["E"] = np.array([[1],[2]])
    res.nodeFields["nodeFields"] = np.arange(4)
    res.nodeFields["nodeFields"].shape = (4,1)

    return "ok"


if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover