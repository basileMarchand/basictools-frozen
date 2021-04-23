# -*- coding: utf-8 -*-

def MeshToPyVista(mesh,TagsAsFields=False):
    from BasicTools.Containers.vtkBridge import MeshToVtk
    import pyvista as pv
    return pv.wrap(MeshToVtk(mesh,TagsAsFields=TagsAsFields))


def PyVistaToMesh(pvmesh,FieldsAsTags=False):
    from BasicTools.Containers.vtkBridge import VtkToMesh
    return VtkToMesh(pvmesh,FieldsAsTags=FieldsAsTags)


def PlotMesh(mesh,**kargs):# pragma: no cover

    from BasicTools.Containers.MeshBase import MeshBase

    if isinstance(mesh,MeshBase):
        pyVistaMesh = MeshToPyVista(mesh)
    else:
        pyVistaMesh = mesh

    pyVistaMesh.plot(**kargs)

def CheckIntegrity(GUI=False):
    from BasicTools.Helpers.Tests import SkipTest
    if SkipTest("PYVISTA_NO_FAIL"): return "ok"

    import BasicTools.Containers.ElementNames as ElementNames
    from BasicTools.Containers.UnstructuredMeshCreationTools import CreateMeshOf
    from BasicTools.Containers.MeshTools import IsClose

    points  = [[0,0,0],[1,0,0],[0,1,0],[0,0,1] ]
    tet = [[0,1,2,3]]
    res = CreateMeshOf(points,tet,elemName = ElementNames.Tetrahedron_4 )

    res.nodeFields = {"x": res.nodes[:,0].flatten(), "Pos":res.nodes}
    res.nodesTags.CreateTag("FirstPoint").AddToTag(0)
    res.elemFields = {"SecondPoint": res.GetElementsOfType(ElementNames.Tetrahedron_4).connectivity[:,1].flatten().astype(float), "conn": res.GetElementsOfType(ElementNames.Tetrahedron_4).connectivity }
    res.GetElementsOfType(ElementNames.Tetrahedron_4).tags.CreateTag("FirstTetrahedron").AddToTag(0)
    sol = MeshToPyVista(res,TagsAsFields= True)
    resII=PyVistaToMesh(sol,FieldsAsTags=True)

    print(res)
    print(resII)

    if not IsClose(res,resII): # pragma: no cover
        raise(Exception("The meshes are not equal"))

    if GUI:
        PlotMesh(res)
        PlotMesh(resII,eye_dome_lighting=True, cpos=[-1, -1, 0.2], color=True)
        
    return 'ok'

if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover
