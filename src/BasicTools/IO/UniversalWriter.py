# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

"""Universal file writer
"""

def InitAllWriters():
    from BasicTools.IO.IOFactory import InitAllWriters as IAW
    IAW()

def WriteMesh(filename,outmesh,binary=None,writer=None):# pragma: no cover

    from BasicTools.IO.IOFactory import CreateWriter
    if writer is None:
        writer = CreateWriter("."+filename.split(".")[-1])
        writer.SetFileName(filename)
        if not (binary is None):
            writer.SetBinary(binary)
        writer.Open()

    PointFields = None
    PointFieldsNames = None
    if hasattr(outmesh,"nodeFields"):
        PointFieldsNames = list(outmesh.nodeFields.keys())
        PointFields = list(outmesh.nodeFields.values())

    CellFields = None
    CellFieldsNames = None
    if hasattr(outmesh,"elemFields"):
        CellFieldsNames = list(outmesh.elemFields.keys())
        CellFields = list(outmesh.elemFields.values())

    writer.Write(outmesh,PointFieldsNames=PointFieldsNames,PointFields=PointFields,CellFieldsNames=CellFieldsNames,CellFields=CellFields )
    writer.Close()

def CheckIntegrity():
    from BasicTools.IO.IOFactory import CreateWriter, InitAllWriters, RegisterWriterClass
    InitAllWriters()

    print(CreateWriter(".geof"))
    print(CreateWriter(".msh"))
    print(CreateWriter(".mesh"))
    print(CreateWriter(".xdmf"))
    try:
        print(CreateWriter("toto"))
    except:
        pass
    else:
        raise (Exception("this must fail " ))

    #normally this class must have the same API as
    #
    from BasicTools.IO.WriterBase import WriterBase as WriterBase

    class MyCustomWriter(WriterBase):
        pass

    RegisterWriterClass(".myDummyExtention",MyCustomWriter,withError=False)

    print(CreateWriter(".myDummyExtention"))

    return "ok"


if __name__ == '__main__':
    print((CheckIntegrity()))# pragma: no cover
