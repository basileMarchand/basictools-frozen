# -*- coding: utf-8 -*-

import numpy as np

import BasicTools.Containers.ElementNames as ElementNames
from BasicTools.Containers.UnstructuredMesh import UnstructuredMesh
from BasicTools.Containers.UnstructuredMeshTools import CreateMeshOfTriangles

#from file vtkCellType.h  of the vtk sources
vtkNameByNumber = {}
vtkNameByNumber[0] = "VTK_EMPTY_CELL"
vtkNameByNumber[1] = "VTK_VERTEX"
vtkNameByNumber[2] = "VTK_POLY_VERTEX"
vtkNameByNumber[3] = "VTK_LINE"
vtkNameByNumber[4] = "VTK_POLY_LINE"
vtkNameByNumber[5] = "VTK_TRIANGLE"
vtkNameByNumber[6] = "VTK_TRIANGLE_STRIP"
vtkNameByNumber[7] = "VTK_POLYGON"
vtkNameByNumber[8] = "VTK_PIXEL"
vtkNameByNumber[9] = "VTK_QUAD"
vtkNameByNumber[10] = "VTK_TETRA"
vtkNameByNumber[11] = "VTK_VOXEL"
vtkNameByNumber[12] = "VTK_HEXAHEDRON"
vtkNameByNumber[13] = "VTK_WEDGE"
vtkNameByNumber[14] = "VTK_PYRAMID"
vtkNameByNumber[15] = "VTK_PENTAGONAL_PRISM"
vtkNameByNumber[16] = "VTK_HEXAGONAL_PRISM"
vtkNameByNumber[21] = "VTK_QUADRATIC_EDGE"
vtkNameByNumber[22] = "VTK_QUADRATIC_TRIANGLE"
vtkNameByNumber[23] = "VTK_QUADRATIC_QUAD"
vtkNameByNumber[36] = "VTK_QUADRATIC_POLYGON"
vtkNameByNumber[24] = "VTK_QUADRATIC_TETRA"
vtkNameByNumber[25] = "VTK_QUADRATIC_HEXAHEDRON"
vtkNameByNumber[26] = "VTK_QUADRATIC_WEDGE"
vtkNameByNumber[27] = "VTK_QUADRATIC_PYRAMID"
vtkNameByNumber[28] = "VTK_BIQUADRATIC_QUAD"
vtkNameByNumber[29] = "VTK_TRIQUADRATIC_HEXAHEDRON"
vtkNameByNumber[30] = "VTK_QUADRATIC_LINEAR_QUAD"
vtkNameByNumber[31] = "VTK_QUADRATIC_LINEAR_WEDGE"
vtkNameByNumber[32] = "VTK_BIQUADRATIC_QUADRATIC_WEDGE"
vtkNameByNumber[33] = "VTK_BIQUADRATIC_QUADRATIC_HEXAHEDRON"
vtkNameByNumber[34] = "VTK_BIQUADRATIC_TRIANGLE"
vtkNameByNumber[35] = "VTK_CUBIC_LINE"
vtkNameByNumber[41] = "VTK_CONVEX_POINT_SET"
vtkNameByNumber[42] = "VTK_POLYHEDRON"
vtkNameByNumber[51] = "VTK_PARAMETRIC_CURVE"
vtkNameByNumber[52] = "VTK_PARAMETRIC_SURFACE"
vtkNameByNumber[53] = "VTK_PARAMETRIC_TRI_SURFACE"
vtkNameByNumber[54] = "VTK_PARAMETRIC_QUAD_SURFACE"
vtkNameByNumber[55] = "VTK_PARAMETRIC_TETRA_REGION"
vtkNameByNumber[56] = "VTK_PARAMETRIC_HEX_REGION"
vtkNameByNumber[60] = "VTK_HIGHER_ORDER_EDGE"
vtkNameByNumber[61] = "VTK_HIGHER_ORDER_TRIANGLE"
vtkNameByNumber[62] = "VTK_HIGHER_ORDER_QUAD"
vtkNameByNumber[63] = "VTK_HIGHER_ORDER_POLYGON"
vtkNameByNumber[64] = "VTK_HIGHER_ORDER_TETRAHEDRON"
vtkNameByNumber[65] = "VTK_HIGHER_ORDER_WEDGE"
vtkNameByNumber[66] = "VTK_HIGHER_ORDER_PYRAMID"
vtkNameByNumber[67] = "VTK_HIGHER_ORDER_HEXAHEDRON"
vtkNameByNumber[68] = "VTK_LAGRANGE_CURVE"
vtkNameByNumber[69] = "VTK_LAGRANGE_TRIANGLE"
vtkNameByNumber[70] = "VTK_LAGRANGE_QUADRILATERAL"
vtkNameByNumber[71] = "VTK_LAGRANGE_TETRAHEDRON"
vtkNameByNumber[72] = "VTK_LAGRANGE_HEXAHEDRON"
vtkNameByNumber[73] = "VTK_LAGRANGE_WEDGE"
vtkNameByNumber[74] = "VTK_LAGRANGE_PYRAMID"

#---------------------------------------------------------------------------
vtkNumberByElementName = {}

vtkNumberByElementName[ElementNames.Point_1] = 1

vtkNumberByElementName[ElementNames.Bar_2] = 3

vtkNumberByElementName[ElementNames.Triangle_3] = 5
vtkNumberByElementName[ElementNames.Quadrangle_4] = 9
vtkNumberByElementName[ElementNames.Tetrahedron_4] = 10

#vtkNumberByElementName[ElementNames.Hexaedron_8] = 11 # voxel
vtkNumberByElementName[ElementNames.Hexaedron_8] = 12
vtkNumberByElementName[ElementNames.Wedge_6] = 13
vtkNumberByElementName[ElementNames.Pyramid_5] = 14

vtkNumberByElementName[ElementNames.Bar_3] = 21
vtkNumberByElementName[ElementNames.Triangle_6] = 22
vtkNumberByElementName[ElementNames.Quadrangle_8] = 23
vtkNumberByElementName[ElementNames.Tetrahedron_10] = 24
vtkNumberByElementName[ElementNames.Hexaedron_20] = 25

elementNameByVtkNumber = {}

for key,vtknumber in vtkNumberByElementName.items():
    elementNameByVtkNumber[vtknumber] = key

elementNameByVtkNumber[11] = ElementNames.Hexaedron_8   #voxel


def MeshToVtk(mesh, vtkobject=None, TagsAsFields=False):


    # From www.vtk;org/wp-content/updloads/2015/04/file-formats.pdf

    try:
        from paraview.vtk import vtkPolyData as vtkPolyData
        from paraview.vtk import vtkUnstructuredGrid as vtkUnstructuredGrid
        from paraview.vtk import vtkPoints
        from paraview.vtk import vtkFloatArray
        from paraview.vtk import vtkIntArray
        from paraview.vtk import vtkIdList
        from paraview.vtk import vtkImageData
    except :
        from vtk import vtkPolyData
        from vtk import vtkUnstructuredGrid
        from vtk import vtkPoints
        from vtk import vtkFloatArray
        from vtk import vtkIntArray
        from vtk import vtkIdList
        from vtk import vtkImageData

    isimagedata = False
    if vtkobject is None:
        if mesh.IsConstantRectilinear():
            output = vtkImageData()
            isimagedata = True
        else:
            usePoly = True
            for  elementsname,elementContainer in mesh.elements.items():
                if ElementNames.dimension[elementsname] == 3:
                    usePoly = False
                    break
            if usePoly:
                output = vtkPolyData()
            else:
                output = vtkUnstructuredGrid()

    else:
        output = vtkobject # pragma: no cover

    if isimagedata:
        output.SetDimensions(mesh.GetDimensions())
        output.SetOrigin(mesh.GetOrigin())
        output.SetSpacing(mesh.GetSpacing())
        if (hasattr(mesh,"nodeFields") and len(mesh.nodeFields) ) or (hasattr(mesh,"elemFields") and len(mesh.elemFields) ):
            print("Warning for the moment only the mesh is converted (no field)")
            print('please use the a vtkUnstructuredGrid container to transfert all the fields')
        return output
    else:
        output.Allocate(mesh.GetNumberOfElements())
        ##copy points
        pts = vtkPoints()
        pts.Allocate(mesh.GetNumberOfNodes())
        if mesh.nodes.shape[1] == 3 :
            for p in range(mesh.GetNumberOfNodes()):
                point = mesh.nodes[p,:]
                pts.InsertNextPoint(point[0],point[1],point[2])
        else:
            #2DCase
            for p in range(mesh.GetNumberOfNodes()):
                point = mesh.nodes[p,:]
                pts.InsertNextPoint(point[0],point[1],0.0)

        output.SetPoints(pts)

        for elementsname,elementContainer in mesh.elements.items():
            pointIds = vtkIdList()
            npe = elementContainer.GetNumberOfNodesPerElement()
            pointIds.SetNumberOfIds(npe)
            vtknumber = vtkNumberByElementName[elementsname]
            for e in range(elementContainer.GetNumberOfElements()):
                for i in range(npe):
                    pointIds.SetId(i,elementContainer.connectivity[e,i])
                output.InsertNextCell(vtknumber, pointIds)


    if hasattr(mesh,"nodeFields"):
        for name,data in mesh.nodeFields.items():
            if data is None:
                continue
            #VTK_data = numpy_support.numpy_to_vtk(num_array=np.swapaxes(phi,0,2).ravel(), deep=True, array_type=vtk.VTK_FLOAT)
            #VTK_data.SetName(name)
            if np.size(data) != mesh.GetNumberOfNodes() and np.size(data) != 2*mesh.GetNumberOfNodes() and np.size(data) != 3*mesh.GetNumberOfNodes():
                print("field ("+str(name)+") is not consistent : it has " + str(np.size(data)) +" values and the mesh has " +str(mesh.GetNumberOfNodes())+ " nodes" )
                raise
                continue

            pd = vtkFloatArray()
            pd.SetName(name)
            if len(data.shape) == 1:
                pd.SetNumberOfComponents(1)
            else:
                pd.SetNumberOfComponents(data.shape[1])
            pd.SetNumberOfTuples(mesh.GetNumberOfNodes())

            if len(data.shape) > 1:
              cpt = 0
              for i in range(mesh.GetNumberOfNodes()):
                 for j in range(data.shape[1]):
                    pd.SetValue(cpt, data[i,j])
                    cpt +=1
              output.GetPointData().AddArray(pd)
            else:
              cpt = 0
              for i in range(mesh.GetNumberOfNodes()):
                    pd.SetValue(cpt, data[i])
                    cpt +=1
              output.GetPointData().AddArray(pd)

    if TagsAsFields:
        for tag in mesh.nodesTags:
            pd = vtkIntArray()
            pd.SetName(tag.name)
            pd.SetNumberOfComponents(1)
            pd.SetNumberOfTuples(mesh.GetNumberOfNodes())
            pd.FillComponent(0,0);

            for i in tag.GetIds():
                pd.SetValue(i,1 )
            output.GetPointData().AddArray(pd)




    if hasattr(mesh,"elemFields"):
        for name,data in mesh.elemFields.items():
            if data is None:
                continue

            if np.size(data) != mesh.GetNumberOfElements() and np.size(data) != 3*mesh.GetNumberOfElements():
                print("field ("+str(name)+") is not consistent : it has " + str(np.size(data)) +" values and the mesh has " +str(mesh.GetNumberOfElements())+ " elements" )
                continue
            pd = vtkFloatArray()
            pd.SetName(name)

            if len(data.shape) == 1:
                pd.SetNumberOfComponents(1)
            else:
                pd.SetNumberOfComponents(data.shape[1])

            pd.SetNumberOfTuples(mesh.GetNumberOfElements())

            if len(data.shape) > 1:
              cpt = 0
              for i in range(mesh.GetNumberOfElements()):
                 for j in range(data.shape[1]):
                    pd.SetValue(cpt, data[i,j])
                    cpt +=1
            else:
              cpt = 0
              for i in range(mesh.GetNumberOfElements()):
                    pd.SetValue(cpt, data[i])
                    cpt +=1

            output.GetCellData().AddArray(pd)
            #output.GetCellData().SetScalars(pd)

    if TagsAsFields:
        elementTags = mesh.GetNamesOfElemTags()
        for tagname in elementTags:
            ids = mesh.GetElementsInTag(tagname)
            pd = vtkIntArray()
            pd.SetName(tagname)
            pd.SetNumberOfComponents(1)
            pd.SetNumberOfTuples(mesh.GetNumberOfElements())
            pd.FillComponent(0,0);
            for i in ids:
                pd.SetValue(i,1 )
            output.GetCellData().AddArray(pd)

    return output

def VtkToMesh(vtkmesh, meshobject=None, TagsAsFields=False):

    if meshobject is None:
        out = UnstructuredMesh()
    else:

        out = meshobject

    from vtk.util import numpy_support
    data = vtkmesh.GetPoints().GetData()
    out.nodes = numpy_support.vtk_to_numpy(data)

    import numpy as np

    out.originalIDNodes = np.arange(out.GetNumberOfNodes())
    nc = vtkmesh.GetNumberOfCells()



    for i in range(nc):
        cell= vtkmesh.GetCell(i)
        ct = cell.GetCellType()
        et = elementNameByVtkNumber[ct]
        nps = cell.GetNumberOfPoints()
        #polyline case
        # we have to be careful because we potentialy change the number of
        # elements in the mesh if we have polylines
        if ct ==4:
            for j in range(nps-1):
                out.GetElementsOfType(et).AddNewElement([cell.GetPointId(j),cell.GetPointId(j+1) ] ,i)
        else:
            out.GetElementsOfType(et).AddNewElement([cell.GetPointId(j) for j in range(nps)] ,i)
    out.PrepareForOutput()

    if vtkmesh.GetPointData().GetNumberOfArrays():
        for f in range(vtkmesh.GetPointData().GetNumberOfArrays()):
            data =  vtkmesh.GetPointData().GetArray(f)
            name = data.GetName()
            nbcomponnents =  data.GetNumberOfComponents()
            nbtuples  = data.GetNumberOfTuples()

            field = np.empty((nbtuples,nbcomponnents),dtype=float)
            cpt =0
            for i in range(nbtuples):
                 for j in range(nbcomponnents):
                    field[i,j] = data.GetValue(cpt)
                    cpt +=1
            out.nodeFields[name] = field

    EOIds = out.GetElementsOriginalIDs()
    EOIds = np.argsort(EOIds)
    if vtkmesh.GetCellData().GetNumberOfArrays():
        for f in range(vtkmesh.GetCellData().GetNumberOfArrays()):
            data =  vtkmesh.GetCellData().GetArray(f)
            if data is None:
                continue
            name = data.GetName()
            nbcomponnents =  data.GetNumberOfComponents()
            nbtuples  = data.GetNumberOfTuples()

            field = np.empty((nbtuples,nbcomponnents),dtype=float)
            cpt =0
            for i in range(nbtuples):
                 for j in range(nbcomponnents):
                    field[EOIds[i],j] = data.GetValue(cpt)
                    cpt +=1
            out.elemFields[name] = field

    return out

def VtkToMeshMultiblock(vtkObject,OP=VtkToMesh)
    if input.IsA("vtkMultiBlockDataSet"):
        res = list()
        nb = input.GetNumberOfBlock()
        for i in range(nb):
            block = input.GetBlock(i)
            res.append(VtkToMeshMultiblock(block,OP=OP))
    else:
      return OP(input)



def CheckIntegrity_VtkToMesh(GUI=False):
    res = CreateMeshOfTriangles([[0,0,0],[1,0,0],[0,1,0],[0,0,1] ], [[0,1,2],[0,2,3]])
    res.nodeFields = {"x": res.nodes[:,0], "Pos":res.nodes}
    res.nodesTags.CreateTag("FirstPoint").AddToTag(0)
    res.elemFields = {"firstPoint": res.GetElementsOfType(ElementNames.Triangle_3).connectivity[:,0], "conn": res.GetElementsOfType(ElementNames.Triangle_3).connectivity }
    res.GetElementsOfType(ElementNames.Triangle_3).tags.CreateTag("FirstTriangle").AddToTag(0)
    sol = MeshToVtk(res,TagsAsFields= True)

    print("CheckIntegrity_VtkToMesh :")
    print(res)
    print(VtkToMesh(sol))
    return 'ok'

def CheckIntegrity_MeshToVtk(GUI=False):
    res = CreateMeshOfTriangles([[0,0,0],[1,0,0],[0,1,0],[0,0,1] ], [[0,1,2],[0,2,3]])
    res.nodeFields = {"x": res.nodes[:,0], "Pos":res.nodes}
    res.nodesTags.CreateTag("FirstPoint").AddToTag(0)
    res.elemFields = {"firstPoint": res.GetElementsOfType(ElementNames.Triangle_3).connectivity[:,0], "conn": res.GetElementsOfType(ElementNames.Triangle_3).connectivity }
    res.GetElementsOfType(ElementNames.Triangle_3).tags.CreateTag("FirstTriangle").AddToTag(0)
    sol = MeshToVtk(res,TagsAsFields= True)

    res = CreateMeshOfTriangles([[0,0],[1,0],[0,1],[1,1] ], [[0,1,2],[2,1,3]])

    sol = MeshToVtk(res )
    print(sol)
    return "OK"


def CheckIntegrity(GUI=False):
    CheckIntegrity_MeshToVtk(GUI)
    CheckIntegrity_VtkToMesh(GUI)
    return 'ok'

if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover
