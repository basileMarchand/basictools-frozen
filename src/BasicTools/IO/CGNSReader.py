# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

""" CGNS mesh file reader

"""
import os

from BasicTools.Bridges.CGNSBridge import CGNSToMesh
from BasicTools.Helpers.BaseOutputObject import BaseOutputObject

try:
    from CGNS.MAP import load
    cgnsLoaded = True
except:
    cgnsLoaded = False


def ReadCGNS(fileName, time  =None, baseNumberOrName=0, zoneNumberOrName=0):
    """Read a CGNS File from disk

    Parameters
    ----------
    fileName : str
        filename of the file to be read
    time : float, optional
        not coded yet, by default None
    baseNumberOrName : int or str, optional
        name of the base to use, by default 0 (first)
    zoneNumberOrName : int or str, optional
        name of the zone to be read, by default 0 (first)

    Returns
    -------
    UnstructuredMesh
        a BasicTools UnstructuredMesh
    """
    reader = CGNSReader()
    reader.SetFileName(fileName)
    reader.baseNumberOrName = baseNumberOrName
    reader.zoneNumberOrName = zoneNumberOrName
    reader.SetTimeToRead(time)
    res = reader.Read(fileName=fileName)
    return res


class CGNSReader(BaseOutputObject):
    def __init__(self):
        super().__init__()
        self.fileName = None
        self.fieldName = None
        self.baseNumberOrName = None
        self.zoneNumberOrName = None
        self.timeToRead = -1

        self.encoding = None
        self.canHandleTemporal = False

    def SetFileName(self,fileName):

        self.fileName = fileName
        if fileName is None :
            self.__path = None
        else:
            self.filePath = os.path.abspath(os.path.dirname(fileName))+os.sep

    def SetTimeToRead(self, time=None):

        if time is None:
            self.timeToRead = 0.
        else:
            raise Exception("not coded yet")
            self.timeToRead = time


    def Read(self, fileName=None, time=None, baseNumberOrName=0, zoneNumberOrName=0):

        if fileName is not None:
            self.SetFileName(fileName)

        self.SetTimeToRead(time)

        pyTree =  load(self.fileName)[0]
        res = CGNSToMesh(pyTree, baseNumberOrName, zoneNumberOrName)

        res.PrepareForOutput()

        return res


if cgnsLoaded:
    from BasicTools.IO.IOFactory import RegisterReaderClass
    RegisterReaderClass(".cgns", CGNSReader)


def CheckIntegrity(GUI=False):

    import BasicTools.IO.CGNSWriter as CW
    CW.CheckIntegrity()

    from BasicTools.Helpers.Tests import TestTempDir
    tempdir = TestTempDir.GetTempPath()

    mesh = ReadCGNS(fileName = tempdir+os.sep+"toto.cgns", baseNumberOrName = 0, zoneNumberOrName = 0)

    print("Read mesh from cgns:", mesh)

    return 'ok'

if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover
