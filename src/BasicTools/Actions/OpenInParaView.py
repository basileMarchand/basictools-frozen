# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#
                       

from BasicTools.Helpers.Tests import TestTempDir
from BasicTools.IO.PathControler import PathControler as PC

paraviewExec= "paraview"
def OpenInParaView( mesh=None,filename=None, run=True):
    from subprocess import Popen
    from BasicTools.Helpers.which import which

    if which(paraviewExec) is None:
        print("App " + str( paraviewExec) + " Not Found" )

    if filename is None and mesh is None and run:
        Popen([paraviewExec])
        return

    if filename is None:
        from BasicTools.Helpers.Tests import GetUniqueTempFile
        (fd,filename) = GetUniqueTempFile(suffix=".xmf",prefix="ExportedDataBasictools_")
    else:
        filename = PC.GetFullFilenameOnTempDirectory(filename)
        #tempdir = TestTempDir.GetTempPath()
        #filename = tempdir + filename

    if mesh is not None:
        from BasicTools.IO.XdmfWriter import WriteMeshToXdmf
        PointFieldsNames = list(mesh.nodeFields.keys())
        PointFields =  list(mesh.nodeFields.values())
        CellFieldsNames = list(mesh.elemFields.keys())
        CellFields = list(mesh.elemFields.values())

        WriteMeshToXdmf(filename,mesh,PointFieldsNames=PointFieldsNames,PointFields=PointFields,CellFieldsNames=CellFieldsNames,CellFields=CellFields  )

    if run:# pragma: no cover
         Popen([paraviewExec,filename])

def CheckIntegrity(GUI=False):
    from BasicTools.Containers.UnstructuredMeshTools import CreateCube
    mesh = CreateCube(dimensions=[20,21,22],spacing=[2.,2.,2.],ofTetras=True)

    from BasicTools.Actions.OpenInParaView import OpenInParaView
    OpenInParaView(mesh,run=GUI)
    OpenInParaView(mesh,filename="CheckIntegrity_OpenInParaview.xmf",run=GUI)

    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity(GUI=True))
    print("Done")