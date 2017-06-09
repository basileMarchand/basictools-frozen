# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-

import BasicTools.FE.ElementNames as ElementNames


def WriteMesh(filename,out,binary=False):# pragma: no cover
    extention = filename.split(".")[-1].lower()


    if extention ==  "geof":
        from BasicTools.IO.GeofWriter import WriteMeshToGeof
        return WriteMeshToGeof(filename,out)
    elif extention ==  "gmsh":
        from  BasicTools.IO.GmshWriter import WriteMeshToGmsh
        return WriteMeshToGmsh(filename,out)
    elif extention ==  "mesh":
        import BasicTools.IO.MeshWriter as WriteMesh
        return WriteMesh(filename,out=out,binary=binary)


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

    return "ok"


if __name__ == '__main__':
    print((CheckIntegrity()))# pragma: no cover