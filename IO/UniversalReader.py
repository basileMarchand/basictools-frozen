# -*- coding: utf-8 -*-

import OTTools.FE.ElementNames as ElementNames


vtknumbers = {}
vtknumbers[ElementNames.Bar_2] = 3

vtknumbers[ElementNames.Triangle_3] = 5
vtknumbers[ElementNames.Quadrangle_4] = 9
vtknumbers[ElementNames.Hexaedron_8] = 12
vtknumbers[ElementNames.Bar_3] = 21
vtknumbers[ElementNames.Triangle_6] = 22
vtknumbers[ElementNames.Tetrahedron_4] = 10
vtknumbers[ElementNames.Tetrahedron_10] = 24


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

    return MeshToVtk(mesh, vtkobject)


def MeshToVtk(mesh, vtkobject=None):

    try:
        from paraview import vtk
    except :
        import vtk

    if vtkobject is not None:
        output = vtkobject
    else:
        output = vtk.vtkUnstructuredGrid()


    output.Allocate(mesh.GetNumberOfElements())
    ##copy points
    pts = vtk.vtkPoints()
    pts.Allocate(mesh.GetNumberOfNodes())
    if mesh.nodes.shape[1] == 3 :
        for p in xrange(mesh.GetNumberOfNodes()):
            point = mesh.nodes[p,:]
            pts.InsertNextPoint(point[0],point[1],point[2])
    else:
        #2DCase
        for p in xrange(mesh.GetNumberOfNodes()):
            point = mesh.nodes[p,:]
            pts.InsertNextPoint(point[0],point[1],0.0)

    output.SetPoints(pts)

    if hasattr(mesh,"nodeFields"):
        for name,data in mesh.nodeFields.iteritems():
            pd = vtk.vtkFloatArray()
            pd.SetName(name)
            pd.SetNumberOfComponents(data.shape[1])
            pd.SetNumberOfTuples(mesh.GetNumberOfNodes())
            cpt = 0
            for i in xrange(mesh.GetNumberOfNodes()):
                for j in xrange(data.shape[1]):
                    pd.SetValue(cpt, data[i,j])
                    cpt +=1

            output.GetPointData().AddArray(pd)

    for elementsname,elementContainer in mesh.elements.iteritems():
        pointIds = vtk.vtkIdList()
        npe = elementContainer.GetNumberOfNodesPerElement()
        pointIds.SetNumberOfIds(npe)
        vtknumber = vtknumbers[elementsname]
        for e in xrange(elementContainer.GetNumberOfElements()):
            for i in xrange(npe):
                pointIds.SetId(i,elementContainer.connectivity[e,i])
            output.InsertNextCell(vtknumber, pointIds)

    if hasattr(mesh,"elemFields"):
        for name,data in mesh.elemFields.iteritems():
            pd = vtk.vtkFloatArray()
            pd.SetName(name)
            pd.SetNumberOfComponents(data.shape[1])
            pd.SetNumberOfTuples(mesh.GetNumberOfElements())
            cpt = 0
            for i in xrange(mesh.GetNumberOfElements()):
                for j in xrange(data.shape[1]):
                    pd.SetValue(cpt, data[i,j])
                    cpt +=1

            output.GetCellData().AddArray(pd)
    return output


def CheckIntegrity():
    from OTTools.FE.UnstructuredMeshTools import CreateMeshOfTriangles
    import numpy as np

    res = CreateMeshOfTriangles([[0,0,0],[1,0,0],[0,1,0],[0,0,1] ], [[0,1,2],[0,2,3]])
    res.elemFields["E"] = np.array([[1],[2]])
    res.nodeFields["nodeFields"] = np.arange(4)
    res.nodeFields["nodeFields"].shape = (4,1)
    sol = MeshToVtk(res)
    print(sol)
    return "ok"


if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover