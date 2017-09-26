# -*- coding: utf-8 -*-

__author__ = "Felipe Bordeu"
import BasicTools.FE.ElementNames as ElementNames


def ReadMesh(filename,out=None):# pragma: no cover
    import os.path

    dirname = os.path.dirname(filename)
    basename,extention = os.path.splitext(os.path.basename(filename))

    if extention ==  ".asc":
        import BasicTools.IO.AscReader as AscReader
        return AscReader.ReadAsc(filename)
    elif extention ==  ".geof":
        import BasicTools.IO.GeofReader as GeofReader
        return GeofReader.ReadGeof(filename)
    elif extention ==  ".msh":
        import BasicTools.IO.GmshReader as GmshReader
        return GmshReader.ReadGmsh(filename,out=out)
    elif extention ==  ".inp":
        import BasicTools.IO.InpReader as ImpReader
        return ImpReader.ReadInp(filename)
    elif extention ==  ".mesh":
        import BasicTools.IO.MeshReader as ReadMesh
        return ReadMesh.ReadMesh(fileName=filename)
    elif extention ==  ".meshb":
        import BasicTools.IO.MeshReader as MeshReader
        return MeshReader.ReadMesh(fileName=filename)
    elif extention ==  ".gcode":
        import BasicTools.IO.GReader as GReader
        return GReader.ReadGCode(fileName=filename)
    elif extention ==  ".fem":
        import BasicTools.IO.FemReader as FemReader
        return FemReader.ReadFem(fileName=filename)
    elif extention ==  ".solb" or extention ==  ".sol":
        import BasicTools.IO.MeshReader as MeshReader

        if extention[-1] == "b":
            f = os.path.join(dirname,basename+".meshb")
        else:
            f = os.path.join(dirname,basename+".mesh")

        # we check if the file exist, if not we try the other type
        if not os.path.isfile(f):
            if extention[-1] == "b":
                f = os.path.join(dirname,basename+".mesh")
            else:
                f = os.path.join(dirname,basename+".meshb")

        if not os.path.isfile(f):
            raise Exception("unable to find a mesh file")
        reader = MeshReader.MeshReader()
        reader.SetFileName(fileName=f)
        reader.Read()
        mesh = reader.output
        fields = reader.ReadExtraFields(filename);
        mesh.nodeFields = {k:v for k,v in fields.items() if k.find("SolAtVertices") != -1  }
        if 'SolAtTetrahedra0' in fields:

            if mesh.GetElementsOfType(ElementNames.Tetrahedron_4).GetNumberOfElements() == mesh.GetNumberOfElements():
                mesh.elemFields = {k:v for k,v in fields.items() if k.find("SolAtTetrahedra") != -1  }
        return mesh
    elif extention ==  ".ut" or extention ==  ".utp":
        from BasicTools.IO.UtReader import UtReader
        reader = UtReader()
        reader.SetFileName(fileName=filename)
        reader.ReadMetaData()
        import BasicTools.IO.GeofReader as GeofReader
        mesh = GeofReader.ReadGeof(reader.meshfile)

        #mesh.nodeFields = {}
        nodesfields = reader.node
        nodesfields.extend(reader.integ)
        for nf in nodesfields:
            mesh.nodeFields[nf] = reader.Read(fieldname=nf,time=-1)
        return mesh
    else:
        raise Exception ("Unkown file extention : " + str(extention))


## to use this function add this lines to the
## programmmable source in paraview, and change the output type to UnstructuredGrid
##
#
#UniversalReader.ReadMeshAnPopulateVtkObject(filename,self.GetOutput())
#from BasicTools.IO.UniversalReader import ReadMeshAndPopulateVtkObject as ReadMeshAndPopulateVtkObject
#filename = "here you put your filename"
#ReadMeshAndPopulateVtkObject(filename,self.GetOutput(),TagsAsFields=True/False )

def ReadMeshAndPopulateVtkObject(filename, vtkobject= None,TagsAsFields=False):# pragma: no cover
    mesh = ReadMesh(filename)
    from BasicTools.FE.UnstructuredMeshTools import MeshToVtk
    return MeshToVtk(mesh, vtkobject,TagsAsFields=TagsAsFields)

def CheckIntegrity():

    return "ok"


if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover