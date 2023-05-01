# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

"""Universal file reader
"""

def InitAllReaders():
    from BasicTools.IO.IOFactory import InitAllReaders as IAR
    IAR()

def ReadMesh(filename,out=None,timeToRead=-1):# pragma: no cover
    import os.path

    dirname = os.path.dirname(filename)
    basename,extention = os.path.splitext(os.path.basename(filename))

    from BasicTools.IO.IOFactory import CreateReader


    reader = CreateReader("."+filename.split(".")[-1])
    reader.SetFileName(filename)
    if reader.canHandleTemporal :
        reader.ReadMetaData()
        reader.SetTimeToRead(timeToRead)
        if timeToRead == -1:
            print("Reading last available time step")
        else:
            print("Reading Time")
            print(timeToRead)

    return reader.Read()

def ReadMeshAndPopulateVtkObject(filename, vtkobject= None,TagsAsFields=False):# pragma: no cover
    mesh = ReadMesh(filename)
    from BasicTools.Bridges.vtkBridge import MeshToVtk
    return MeshToVtk(mesh, vtkobject,TagsAsFields=TagsAsFields)

def CheckIntegrity():
    return "ok"


if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover
