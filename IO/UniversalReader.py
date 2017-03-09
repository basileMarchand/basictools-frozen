# -*- coding: utf-8 -*-

import OTTools.FE.ElementNames as ElementNames




def ReadMesh(filename,out=None):# pragma: no cover
    extention = filename.split(".")[-1].lower()
    if extention ==  "asc":
        import OTTools.IO.AscReader as AscReader
        return AscReader.ReadAsc(filename)
    elif extention ==  "geof":
        import OTTools.IO.GeofReader as GeofReader
        return GeofReader.ReadGeof(filename)
    elif extention ==  "msh":
        import OTTools.IO.GmshReader as GmshReader
        return GmshReader.ReadGmsh(filename,out=out)
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
    elif extention ==  "fem":
        import OTTools.IO.FemReader as FemReader

        return FemReader.ReadFem(fileName=filename)

    elif extention ==  "solb" or extention ==  "sol":
        import OTTools.IO.MeshReader as MeshReader
        if extention[-1] == "b":
            f = filename.split(".")[0] + ".meshb"
        else:
            f = filename.split(".")[0] + ".mesh"

        # we check if the file exist, if not we try the other type
        import os.path
        if not os.path.isfile(f):
            if extention[-1] == "b":
                f = filename.split(".")[0] + ".mesh"
            else:
                f = filename.split(".")[0] + ".meshb"


        reader = MeshReader.MeshReader()
        reader.SetFileName(fileName=f)
        reader.Read()
        mesh = reader.output
        fields = reader.ReadExtraField(filename);
        mesh.nodeFields = {k:v for k,v in fields.iteritems() if k.find("SolAtVertices") != -1  }
        if fields.has_key('SolAtTetrahedra0'):

            if mesh.GetElementsOfType(ElementNames.Tetrahedron_4).GetNumberOfElements() == mesh.GetNumberOfElements():
                mesh.elemFields = {k:v for k,v in fields.iteritems() if k.find("SolAtTetrahedra") != -1  }
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
#ReadMeshAndPopulateVtkObject(filename,self.GetOutput(),TagsAsFields=True/False )

def ReadMeshAndPopulateVtkObject(filename, vtkobject= None,TagsAsFields=False):# pragma: no cover
    mesh = ReadMesh(filename)
    from OTTools.FE.UnstructuredMeshTools import MeshToVtk
    return MeshToVtk(mesh, vtkobject,TagsAsFields=TagsAsFields)

def CheckIntegrity():

    return "ok"


if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover