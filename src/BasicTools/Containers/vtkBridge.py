# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
                       

import numpy as np

import BasicTools.Containers.ElementNames as ElementNames
from BasicTools.Containers.UnstructuredMesh import UnstructuredMesh
from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshOfTriangles

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
vtkNumberByElementName[ElementNames.Quadrangle_9] = 28
elementNameByVtkNumber = {}

for key,vtknumber in vtkNumberByElementName.items():
    elementNameByVtkNumber[vtknumber] = key

elementNameByVtkNumber[11] = ElementNames.Hexaedron_8   #voxel

def VtkFieldToNumpyField(support,vtkField):
    from vtk.util import numpy_support

    name = vtkField.GetName()

    data = numpy_support.vtk_to_numpy(vtkField)
    if support.IsConstantRectilinear():
        dims = list(support.GetDimensions())[::-1]
        dims.append(data.shape[1])
        data.shape = tuple(dims)
        data = np.swapaxes(data,0,2)

    return (name,data)

def NumpyFieldToVtkField(support,fielddata,fieldname):
    from vtk.util import numpy_support


    isimagedata = support.IsConstantRectilinear()
    outputtype = numpy_support.get_vtk_array_type(fielddata.dtype)

    if len(fielddata.shape) > 1:
      if isimagedata:
          dataView = fielddata.view()
          dims = list(support.GetDimensions())
          dims.append(fielddata.shape[1])
          dataView.shape = tuple(dims)
          dataView.shape = support.GetDimensions()
          VTK_data = numpy_support.numpy_to_vtk(num_array=np.swapaxes(dataView,0,2).ravel(), deep=True, array_type=outputtype)
      else:
          VTK_data = numpy_support.numpy_to_vtk(num_array=fielddata, deep=True, array_type=outputtype)
    else:
      #cpt = 0
      if isimagedata:
          dataView = fielddata.view()
          dataView.shape = support.GetDimensions()
          VTK_data = numpy_support.numpy_to_vtk(num_array=np.swapaxes(dataView,0,2).ravel(), deep=True, array_type=outputtype)
      else:
          VTK_data = numpy_support.numpy_to_vtk(num_array=fielddata, deep=True, array_type=outputtype)
    VTK_data.SetName(fieldname)
    return VTK_data

def ApplyVtkPipeline(mesh,op):
    vtkMesh = MeshToVtk(mesh)
    vtkOuputMesh = op(vtkMesh)
    return VtkToMesh(vtkOuputMesh)

def MeshToVtk(mesh, vtkobject=None, TagsAsFields=False):


    # From www.vtk;org/wp-content/updloads/2015/04/file-formats.pdf

    try:
        from paraview.vtk import vtkPolyData, vtkUnstructuredGrid, vtkPoints,vtkIdList, vtkImageData
    except :
        from vtk import vtkPolyData, vtkUnstructuredGrid, vtkPoints, vtkIdList, vtkImageData

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
        #if (hasattr(mesh,"nodeFields") and len(mesh.nodeFields) ) or (hasattr(mesh,"elemFields") and len(mesh.elemFields) ):
        #    print("Warning for the moment only the mesh is converted (no field)")
        #    print('please use the a vtkUnstructuredGrid container to transfert all the fields')
        #return output
    else:
        output.Allocate(mesh.GetNumberOfElements())
        ##copy points
        pts = vtkPoints()
        pts.Allocate(mesh.GetNumberOfNodes())

        VTK_originalIDNodes = NumpyFieldToVtkField(mesh,mesh.originalIDNodes,"originalIds")
        output.GetPointData().AddArray(VTK_originalIDNodes)


        #nodeOriginalIds = vtkIntArray()
        #nodeOriginalIds.SetName("originalIds")
        #nodeOriginalIds.SetNumberOfComponents(1)
        #nodeOriginalIds.SetNumberOfTuples(mesh.GetNumberOfNodes())
        if mesh.nodes.shape[1] == 3 :
            for p in range(mesh.GetNumberOfNodes()):
                point = mesh.nodes[p,:]
                pts.InsertNextPoint(point[0],point[1],point[2])
                #nodeOriginalIds.SetValue(p, mesh.originalIDNodes[p])
        else:
            #2DCase
            for p in range(mesh.GetNumberOfNodes()):
                point = mesh.nodes[p,:]
                pts.InsertNextPoint(point[0],point[1],0.0)
                #nodeOriginalIds.SetValue(p, mesh.originalIDNodes[p])

        #output.GetPointData().AddArray(nodeOriginalIds)
        output.SetPoints(pts)

        VTK_originalIDsEl = NumpyFieldToVtkField(mesh,mesh.GetElementsOriginalIDs(),"originalIds")
        output.GetCellData().AddArray(VTK_originalIDsEl)

        #elemOriginalIds = vtkIntArray()
        #elemOriginalIds.SetName("originalIds")
        #elemOriginalIds.SetNumberOfComponents(1)
        #elemOriginalIds.SetNumberOfTuples(mesh.GetNumberOfElements())
        cpt = 0
        for elementsname,elementContainer in mesh.elements.items():
            pointIds = vtkIdList()
            npe = elementContainer.GetNumberOfNodesPerElement()
            pointIds.SetNumberOfIds(npe)
            vtknumber = vtkNumberByElementName[elementsname]
            for e in range(elementContainer.GetNumberOfElements()):
                for i in range(npe):
                    pointIds.SetId(i,elementContainer.connectivity[e,i])
                output.InsertNextCell(vtknumber, pointIds)
                #elemOriginalIds.SetValue(cpt, elementContainer.originalIds[e])
                cpt += 1
        #output.GetCellData().AddArray(elemOriginalIds)

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

            #print("name : ", name )
            VTK_data = NumpyFieldToVtkField(mesh,data,name)
            output.GetPointData().AddArray(VTK_data)
            continue

    if TagsAsFields:
        tagMask = np.empty(mesh.GetNumberOfNodes(),int)

        for tag in mesh.nodesTags:
            tag.GetIdsAsMask(output=tagMask)
            VTK_data = NumpyFieldToVtkField(mesh,tagMask,tag.name)
            output.GetPointData().AddArray(VTK_data)
            continue

    if hasattr(mesh,"elemFields"):
        for name,data in mesh.elemFields.items():

            if data is None:
                continue
            if np.size(data) != mesh.GetNumberOfElements() and np.size(data) != 2*mesh.GetNumberOfElements() and np.size(data) != 3*mesh.GetNumberOfElements():
                print("field ("+str(name)+") is not consistent : it has " + str(np.size(data)) +" values and the mesh has " +str(mesh.GetNumberOfNodes())+ " nodes" )
                raise
                continue

            #print("name : ", name )
            VTK_data = NumpyFieldToVtkField(mesh,data,name)
            output.GetCellData().AddArray(VTK_data)
            continue

    if TagsAsFields:
        elementTags = mesh.GetNamesOfElemTags()
        #print(elementTags)
        for tagname in elementTags:
            ids = mesh.GetElementsInTag(tagname)
            tagMask = np.zeros(mesh.GetNumberOfElements(),int )
            tagMask[ids] = True
            VTK_data = NumpyFieldToVtkField(mesh,tagMask,tagname)
            output.GetCellData().AddArray(VTK_data)
            continue

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
        if ct == 4:
            for j in range(nps-1):
                out.GetElementsOfType(et).AddNewElement([cell.GetPointId(j),cell.GetPointId(j+1) ] ,i)
        elif ct ==  11:
            # 11 is a voxel and the numbering is not the same as the hexahedron
            #https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
            original_coonectivity = [cell.GetPointId(j) for j in range(nps)]
            connectivity = original_coonectivity[[0,1,3,24,5,7,6]]
            out.GetElementsOfType(et).AddNewElement(connectivity  ,i)
        else:
            out.GetElementsOfType(et).AddNewElement([cell.GetPointId(j) for j in range(nps)] ,i)
    out.PrepareForOutput()

    if vtkmesh.GetPointData().GetNumberOfArrays():
        for f in range(vtkmesh.GetPointData().GetNumberOfArrays()):
            data =  vtkmesh.GetPointData().GetArray(f)
            (name,field) = VtkFieldToNumpyField(out,data)
            if name == "originalIds":
                out.originalIds = field
            else:
                if len(field.shape) == 1 and field.dtype ==  int and min(field)>=0 and max(field) <= 1:
                    out.nodesTags.CreateTag(name).SetIds(np.where(field)[0])
                else:
                    out.nodeFields[name] = field

    EOIds = out.GetElementsOriginalIDs()
    EOIds = np.argsort(EOIds)
    if vtkmesh.GetCellData().GetNumberOfArrays():
        for f in range(vtkmesh.GetCellData().GetNumberOfArrays()):
            data =  vtkmesh.GetCellData().GetArray(f)
            if data is None:
                continue

            data =  vtkmesh.GetCellData().GetArray(f)
            (name,field) = VtkFieldToNumpyField(out,data)
            Elfield = np.empty(field.shape,dtype=float)
            if len(field.shape) > 1:
                Elfield[EOIds,:] = field[range(field.shape[0]),:]
            else:
                Elfield[EOIds] = field[:]
            if name == "originalIds":
                out.SetElementsOriginalIDs(Elfield)
            else:
                if len(field.shape) == 1 and field.dtype ==  int and min(field)>=0 and max(field) <= 1 :
                    cpt = 0
                    for elname,data in out.elements.items():
                        nn = data.GetNumberOfElements()
                        data.tags.CreateTag(name).SetIds(np.where(field[cpt:cpt+nn])[0])
                        cpt += nn

                else:
                    out.elemFields[name] = field
            continue

    return out

def VtkToMeshMultiblock(vtkObject,OP=VtkToMesh):
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
    res.nodeFields = {"x": res.nodes[:,0].flatten(), "Pos":res.nodes}
    res.nodesTags.CreateTag("FirstPoint").AddToTag(0)
    res.elemFields = {"SecondPoint": res.GetElementsOfType(ElementNames.Triangle_3).connectivity[:,1].flatten(), "conn": res.GetElementsOfType(ElementNames.Triangle_3).connectivity }
    res.GetElementsOfType(ElementNames.Triangle_3).tags.CreateTag("FirstTriangle").AddToTag(0)
    sol = MeshToVtk(res,TagsAsFields= True)

    print("CheckIntegrity_VtkToMesh :")
    print(res)
    resII=VtkToMesh(sol)
    print(resII)
    from BasicTools.Containers.MeshTools import IsClose
    if not IsClose(res,resII):
        raise(Exception("The meshes are not equal"))
    return 'ok'

def CheckIntegrity_MeshToVtk(GUI=False):
    res = CreateMeshOfTriangles([[0,0,0],[1,0,0],[0,1,0],[0,0,1] ], [[0,1,2],[0,2,3]])
    res.nodeFields = {"x": res.nodes[:,0].flatten(), "Pos":res.nodes}
    res.nodesTags.CreateTag("FirstPoint").AddToTag(0)
    res.elemFields = {"SecondPoint": res.GetElementsOfType(ElementNames.Triangle_3).connectivity[:,1].flatten(), "conn": res.GetElementsOfType(ElementNames.Triangle_3).connectivity }
    res.GetElementsOfType(ElementNames.Triangle_3).tags.CreateTag("FirstTriangle").AddToTag(0)
    sol = MeshToVtk(res,TagsAsFields= True)
    print(sol)
    resII = VtkToMesh(sol)
    print(res)
    print(resII)
    from BasicTools.Containers.MeshTools import IsClose
    if not IsClose(res,resII):
        raise(Exception("The meshes are not equal"))

    ## test a 2D mesh
    res = CreateMeshOfTriangles([[0,0],[1,0],[0,1],[1,1] ], [[0,1,2],[2,1,3]])
    sol = MeshToVtk(res )


    return "OK"


def CheckIntegrity(GUI=False):
    CheckIntegrity_MeshToVtk(GUI)
    CheckIntegrity_VtkToMesh(GUI)
    return 'ok'

if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover
