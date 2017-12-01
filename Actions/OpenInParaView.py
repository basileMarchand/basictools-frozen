# -*- coding: utf-8 -*-

from BasicTools.Helpers.Tests import TestTempDir
from BasicTools.IO.PathControler import PathControler as PC

paraviewExec= "paraview"
def OpenInParaView( mesh=None,filename=None ):

    if filename is None:
        from BasicTools.Helpers.Tests import GetUniqueTempFile
        (fd,filename) = GetUniqueTempFile(suffix=".xmf",prefix="ExportedDataBasictools_")
    else:
        filename = PC.GetFullFilenameOnTempDirectory(filename)
        #tempdir = TestTempDir.GetTempPath()
        #filename = tempdir + filename

    if mesh is not None:
        from BasicTools.IO.XdmfWriter import WriteMeshToXdmf
        PointFieldsNames = mesh.nodeFields.keys()
        PointFields =  mesh.nodeFields.values()
        CellFieldsNames = mesh.elemFields.keys()
        CellFields = mesh.elemFields.values()

        WriteMeshToXdmf(filename,mesh,PointFieldsNames=PointFieldsNames,PointFields=PointFields,CellFieldsNames=CellFieldsNames,CellFields=CellFields  )

    from subprocess import Popen
    Popen([paraviewExec,filename])

def CheckIntegrity():
    pass
