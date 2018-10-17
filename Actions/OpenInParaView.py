# -*- coding: utf-8 -*-

from BasicTools.Helpers.Tests import TestTempDir
from BasicTools.IO.PathControler import PathControler as PC

paraviewExec= "paraview"
def OpenInParaView( mesh=None,filename=None, run=True):

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

    if run:
        from subprocess import Popen
        Popen([paraviewExec,filename])

def CheckIntegrity():
    from BasicTools.Actions.OpenInParaView import OpenInParaView
    #OpenInParaView(mesh)
    pass

