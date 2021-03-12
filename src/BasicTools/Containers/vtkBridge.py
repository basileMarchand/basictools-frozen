# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#


import numpy as np

import BasicTools.Containers.ElementNames as ElementNames
from BasicTools.Containers.UnstructuredMesh import UnstructuredMesh, AllElements
from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshOfTriangles
from BasicTools.TestData import GetTestDataPath

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

vtkNumberByElementName[ElementNames.Pyramid_13] = 27
vtkNumberByElementName[ElementNames.Wedge_15] = 26

elementNameByVtkNumber = {}

for key,vtknumber in vtkNumberByElementName.items():
    elementNameByVtkNumber[vtknumber] = key

elementNameByVtkNumber[4] = ElementNames.Bar_2
elementNameByVtkNumber[11] = ElementNames.Hexaedron_8   #voxel

#if a field is of type [..] and the min max are 0 and 1 then the field is
#converted to a tag. the first type is used to encode tags in vtk
tagsTypes = [np.int8, np.uint8, int]

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

def PlotMesh(mesh):
    import vtk

    from BasicTools.Containers.MeshBase import MeshBase
    if isinstance(mesh,MeshBase):
        vtkMesh = MeshToVtk(mesh)
    else:
        vtkMesh = mesh

    vGF = vtk.vtkGeometryFilter()
    vGF.SetInputData(vtkMesh)
    vGF.Update()

    nbArrays = vtkMesh.GetPointData().GetNumberOfArrays()


    cylinderMapper = vtk.vtkPolyDataMapper()
    if nbArrays > 0:
        out2 = vtk.vtkAssignAttribute()
        out2.SetInputConnection(vGF.GetOutputPort())

        array = vtkMesh.GetPointData().GetArray(0)
        out2.Assign(array.GetName(), "SCALARS", "POINT_DATA")
        cylinderMapper.SetInputConnection(out2.GetOutputPort())
    else:
        cylinderMapper.SetInputConnection(vGF.GetOutputPort())

    cylinderActor = vtk.vtkActor()
    cylinderActor.SetMapper(cylinderMapper)
    cylinderActor.RotateX(30.0)
    cylinderActor.RotateY(-45.0)

    ren = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)

    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)

    style = vtk.vtkInteractorStyleTrackballCamera()
    style.SetDefaultRenderer(ren)
    iren.SetInteractorStyle(style)

    ren.AddActor(cylinderActor)
    ren.GradientBackgroundOn();
    ren.SetBackground(0.3176,0.3412,0.4314);
    ren.SetBackground2(0,0,0.1647);
    renWin.SetSize(800, 600)

    buttonWidget = vtk.vtkButtonWidget()
    buttonWidget.SetInteractor(iren)

    class Listener():
        def __init__(self):
            self.fildcpt = 0
            self.minmaxs = dict()

        def processStateChangeEvent(self,obj,ev):

            self.fildcpt += 1
            if vtkMesh.GetPointData().GetNumberOfArrays() == 0:
                return

            nb = self.fildcpt % vtkMesh.GetPointData().GetNumberOfArrays()
            array = vtkMesh.GetPointData().GetArray(nb)
            arrayName = array.GetName()

            out2.Assign(arrayName, "SCALARS", "POINT_DATA")
            if arrayName in self.minmaxs:
               lo, hi =  self.minmaxs[arrayName]
            else:
               lo,hi = array.GetRange()
               self.minmaxs[arrayName] = (lo,hi)

            lut = vtk.vtkColorTransferFunction()
            lut.SetColorSpaceToHSV()
            lut.SetColorSpaceToDiverging()
            lut.AddRGBPoint(lo,0.23137254902000001, 0.298039215686, 0.75294117647100001)
            lut.AddRGBPoint((lo+hi)/2,0.86499999999999999, 0.86499999999999999, 0.86499999999999999)
            lut.AddRGBPoint(hi,0.70588235294099999, 0.015686274509800001, 0.149019607843)
            lut.Build()

            cylinderMapper.SetInterpolateScalarsBeforeMapping(True)
            cylinderMapper.SetScalarRange(lo, hi)
            cylinderMapper.SetUseLookupTableScalarRange(True)
            cylinderMapper.SetLookupTable(lut)
            cylinderMapper.SetScalarModeToUsePointData()
            cylinderMapper.ScalarVisibilityOn()
            cylinderMapper.SelectColorArray(array.GetName())
            print("Plot of field {} min/max, {}/{}".format(arrayName,lo,hi))
            renWin.Render()

            pass
    listener = Listener()
    listener.processStateChangeEvent(None,None)

    buttonWidget.AddObserver( 'StateChangedEvent', listener.processStateChangeEvent )

    buttonRepresentation= vtk.vtkTexturedButtonRepresentation2D()
    buttonRepresentation.SetNumberOfStates(1);
    r = vtk.vtkPNGReader()

    fileName  = GetTestDataPath()+"Next.png"

    r.SetFileName(fileName)
    r.Update()
    image = r.GetOutput()
    buttonRepresentation.SetButtonTexture(0, image)

    buttonRepresentation.SetPlaceFactor(1);
    buttonRepresentation.PlaceWidget([0,50,0,50,0,0]);

    buttonWidget.SetRepresentation(buttonRepresentation);

    buttonWidget.On();

    iren.Initialize()

    axesActor = vtk.vtkAxesActor();
    axes = vtk.vtkOrientationMarkerWidget()
    axes.SetOrientationMarker(axesActor)
    axes.SetInteractor(iren)
    axes.EnabledOn()
    axes.InteractiveOn()
    ren.ResetCamera()
    ren.GetActiveCamera().Zoom(1.5)
    renWin.Render()
    # Start the event loop.
    iren.Start()
    renWin.Finalize()
    iren.TerminateApp()
    del renWin, iren

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
        pts.SetNumberOfPoints(mesh.GetNumberOfNodes())

        VTK_originalIDNodes = NumpyFieldToVtkField(mesh,mesh.originalIDNodes,"originalIds")
        output.GetPointData().AddArray(VTK_originalIDNodes)

        if mesh.nodes.shape[1] == 3 :
            for p in range(mesh.GetNumberOfNodes()):
                point = mesh.nodes[p,:]
                pts.SetPoint(p,point[0],point[1],point[2])
        else:
            for p in range(mesh.GetNumberOfNodes()):
                point = mesh.nodes[p,:]
                pts.SetPoint(p,point[0],point[1],0.0)

        output.SetPoints(pts)

        VTK_originalIDsEl = NumpyFieldToVtkField(mesh,mesh.GetElementsOriginalIDs(),"originalIds")
        output.GetCellData().AddArray(VTK_originalIDsEl)

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
                cpt += 1

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
        tagMask = np.empty(mesh.GetNumberOfNodes(),tagsTypes[0])

        for tag in mesh.nodesTags:
            tag.GetIdsAsMask(output=tagMask)
            VTK_data = NumpyFieldToVtkField(mesh,tagMask,tag.name)
            output.GetPointData().AddArray(VTK_data)
            continue

    if hasattr(mesh,"elemFields"):
        for name,data in mesh.elemFields.items():

            if data is None:
                continue

            if mesh.GetNumberOfElements() == 0:
                continue

            if np.size(data)/mesh.GetNumberOfElements() !=  np.size(data)//mesh.GetNumberOfElements() :
                print("field ("+str(name)+") is not consistent : it has " + str(np.size(data)) +" values and the mesh has " +str(mesh.GetNumberOfElements())+ " elements" )
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
            tagMask = np.zeros(mesh.GetNumberOfElements(),dtype=tagsTypes[0] )
            tagMask[ids] = True
            VTK_data = NumpyFieldToVtkField(mesh,tagMask,tagname)
            output.GetCellData().AddArray(VTK_data)
            continue

    return output

def VtkToMeshOnlyMeta(vtkmesh, FieldsAsTags=False):
    from vtk.util import numpy_support

    class UnstructuredMeshMetaData():
        def __init__(self):
            self.nbnodes = 0
            self.originalIDNodes = False
            self.nodesTags = []
            self.nodeFields= []

            self.nbelements = 0
            self.originalIDElements = False
            self.elemTags = []
            self.elemFields = []

    res = UnstructuredMeshMetaData()
    res.nbnodes = vtkmesh.GetPoints().GetNumberOfPoints()
    res.nbelements = vtkmesh.GetCells().GetNumberOfCells()

    for f in range(vtkmesh.GetCellData().GetNumberOfArrays()):
        data =  vtkmesh.GetCellData().GetAbstractArray(f)

        if data is None:
            continue

        if data.IsNumeric() is None:
            continue

        nptype = numpy_support.get_numpy_array_type(data.GetDataType())
        name = data.GetName()
        if name == "originalIds":
            res.originalIDElements = True
        else:
            rmin,rmax = data.GetRange()
            if FieldsAsTags and  nptype in tagsTypes and rmin>=0 and rmax <= 1:
                res.elemTags.append(name)
            else:
                res.elemFields.append(name)
        continue

    for f in range(vtkmesh.GetPointData().GetNumberOfArrays()):
        data =  vtkmesh.GetPointData().GetAbstractArray(f)

        if data is None:
            continue

        if data.IsNumeric() is None:
            continue

        name = data.GetName()
        nptype = numpy_support.get_numpy_array_type(data.GetDataType())
        if name == "originalIds":
            res.originalIDNodes = True
        else:
            rmin,rmax = data.GetRange()
            if FieldsAsTags and nptype in tagsTypes and rmin>=0 and rmax <= 1 :
                res.nodesTags.append(name)
            else:
                res.nodeFields.append(name)
    return res

def VtkToMesh(vtkmesh, meshobject=None, FieldsAsTags=True):

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
            if nps > 2 :
                print("Warning polyline with more than 2 nodes, elemfield are incompatible after conversion ")
            for j in range(nps-1):
                out.GetElementsOfType(et).AddNewElement([cell.GetPointId(j),cell.GetPointId(j+1) ] ,i)
        elif ct ==  11:
            # 11 is a voxel and the numbering is not the same as the hexahedron
            #https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
            original_coonectivity = np.array([cell.GetPointId(j) for j in range(nps)])
            connectivity = original_coonectivity[[0,1,3,2,4,5,7,6]]

            out.GetElementsOfType(et).AddNewElement(connectivity  ,i)
        else:
            out.GetElementsOfType(et).AddNewElement([cell.GetPointId(j) for j in range(nps)] ,i)
    out.PrepareForOutput()

    if vtkmesh.GetPointData().GetNumberOfArrays():
        for f in range(vtkmesh.GetPointData().GetNumberOfArrays()):
            data =  vtkmesh.GetPointData().GetArray(f)
            (name,field) = VtkFieldToNumpyField(out,data)
            if name == "originalIds":
                out.originalIDNodes = field
            else:
                if FieldsAsTags and len(field.shape) == 1 and field.dtype in tagsTypes and np.min(field)>=0 and np.max(field) <= 1:
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
                if FieldsAsTags and  len(field.shape) == 1 and field.dtype in tagsTypes and np.min(field)>=0 and np.max(field) <= 1 :
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
    if GUI:
        res.nodeFields["Field1"] = np.array([30,20,30,1])
        res.nodeFields["Field2"] = np.array([0,1,0,1])+0.1
        PlotMesh(res)
        sol = MeshToVtk(res )

        PlotMesh(sol)

    return "OK"


def CheckIntegrity(GUI=False):
    CheckIntegrity_MeshToVtk(GUI)
    CheckIntegrity_VtkToMesh(GUI)
    return 'ok'

if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover
