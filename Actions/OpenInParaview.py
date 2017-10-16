# -*- coding: utf-8 -*-

from BasicTools.Helpers.Tests import TestTempDir


def OpenInParaView( mesh,filename=None ):

    if filename is None:
        from BasicTools.Helpers.Tests import GetUniqueTempFile
        (fd,filename) = GetUniqueTempFile(suffix=".xmf",prefix="ExportedDataBasictools_")
    else:
        tempdir = TestTempDir.GetTempPath()
        filename = tempdir + filename


    from BasicTools.IO.XdmfWriter import WriteMeshToXdmf


    PointFieldsNames = mesh.nodeFields.keys()
    PointFields =  mesh.nodeFields.values()
    CellFieldsNames = mesh.elemFields.keys()
    CellFields = mesh.elemFields.values()

    WriteMeshToXdmf(filename,mesh,PointFieldsNames=PointFieldsNames,PointFields=PointFields,CellFieldsNames=CellFieldsNames,CellFields=CellFields  )

    from subprocess import Popen
    Popen(["paraview",filename])

def CheckIntegrity():
    pass
