import os

from BasicTools.Containers.UnstructuredMesh import UnstructuredMesh
from BasicTools.IO.StlWriter import WriteMeshToStl
from BasicTools.Helpers.Tests import TestTempDir, WriteTempFile
from BasicTools.IO.GmshReader import ReadGmsh
from BasicTools.IO.CodeInterface import Interface

gmshExec= "gmsh"
def StlToMesh(stlmesh : UnstructuredMesh) -> UnstructuredMesh:
    """Fill 2D close surface with 3D elements. This functionality uses gmsh whit a external call

    Parameters
    ----------
    stlmesh : UnstructuredMesh
        a close mesh composed only by 2D elements.

    Returns
    -------
    UnstructuredMesh
        a mesh composed with by the 2D exterior mesh filled with 3D elements. The order of the 2D
        element (and nodes) is not necessarily preserved.
    """

    gmshTemplate = """
Merge "{stlFileName}"; // load the surface mesh
Surface Loop(1) = Surface{{:}}; // make a single shell from all surfaces, assuming that the mesh is watertight
Volume(1) = {{1}}; // create a volume from the shell
Mesh.Algorithm3D=10; // test the new 3D algorithm
"""

    inter = Interface(workingDirectory=TestTempDir.GetTempPath())
    inter.inputFileExtension = ".geo"
    inter.inputFilename = "StlToMeshScript"
    inter.tpl = gmshTemplate
    inter.codeCommand = gmshExec

    stlFileName =  "StlToMesh_mesh.stl"
    inter.parameters = {"stlFileName":stlFileName}
    WriteMeshToStl(TestTempDir.GetTempPath() + os.sep+stlFileName,stlmesh)
    inter.WriteFile(0)

    inter.options = ["-3","-f","msh2"]

    print(inter.SingleRunComputationAndReturnOutput(0))

    return ReadGmsh(TestTempDir.GetTempPath() + inter.inputFilename + "0.msh"  )

def CheckIntegrity(GUI=False):
    from BasicTools.Helpers.Tests import SkipTest
    if SkipTest("GMSH_NO_FAIL"):
        return "skip"

    from BasicTools.Helpers.which import which

    if which(gmshExec) is None:
        print("gmsh not found")
        return "skip"

    from BasicTools.IO.StlReader import ReadStl
    from BasicTools.TestData import GetTestDataPath
    stlmesh = ReadStl(GetTestDataPath()+"stlsphere.stl")

    mesh = StlToMesh(stlmesh)
    print(stlmesh)
    print(mesh)
    if GUI:
        from BasicTools.Actions.OpenInParaView import OpenInParaView
        OpenInParaView(mesh)
    return 'ok'

if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover



