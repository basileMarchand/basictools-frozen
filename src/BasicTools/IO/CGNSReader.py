# -*- coding: utf-8 -*-
#
# This file is subject to the terms and conditions defined in
# file 'LICENSE.txt', which is part of this source code package.
#

"""CGNS file reader
"""

import os

import numpy as np

from BasicTools.Bridges.CGNSBridge import CGNSToMesh
from BasicTools.Helpers.BaseOutputObject import BaseOutputObject

def ReadCGNS(fileName, time = None, baseNumberOrName = 0, zoneNumberOrName = 0):
    """Function API for reading a CGNS file

    Parameters
    ----------
    fileName : str
        name of the file to be read
    time : float, optional
        not coded yet, by default None
    baseNumberOrName : int or str, optional
        name of the base to use, by default 0 (first)
    zoneNumberOrName : int or str, optional
        name of the zone to be read, by default 0 (first)

    Returns
    -------
    UnstructuredMesh
        output unstructured mesh object containing reading result
    """
    reader = CGNSReader()
    reader.SetFileName(fileName)
    reader.baseNumberOrName = baseNumberOrName
    reader.zoneNumberOrName = zoneNumberOrName
    reader.SetTimeToRead(time)
    res = reader.Read(fileName=fileName)
    return res


class CGNSReader(BaseOutputObject):
    """CGNS Reader class
    """

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
        """Function to set fileName to read

        Parameters
        ----------
        fileName : str
            name of the file to be read
        """
        self.fileName = fileName
        if fileName is None :
            self.__path = None
        else:
            self.filePath = os.path.abspath(os.path.dirname(fileName))+os.sep

    def SetTimeToRead(self, time=None):
        """Function to set time value to read

        Parameters
        ----------
        time : float, optional
            not coded yet, by default None
        """
        if time is None:
            self.timeToRead = 0.
        else:
            raise Exception("not coded yet")
            self.timeToRead = time


    def Read(self, fileName=None, time=None, baseNumberOrName=0, zoneNumberOrName=0):
        """Function that performs the reading of a CGNS result file

        Parameters
        ----------
        fileName : str
            name of the file to be read
        time : float, optional
            not coded yet, by default None
        baseNumberOrName : int or str, optional
            name of the base to use, by default 0 (first)
        zoneNumberOrName : int or str, optional
            name of the zone to be read, by default 0 (first)

        Returns
        -------
        UnstructuredMesh
            output unstructured mesh object containing reading result
        """
        if fileName is not None:
            self.SetFileName(fileName)

        self.SetTimeToRead(time)


        from h5py import File
        h5file =  File(self.fileName,"r")

        def ConvertData(node):
            res = [node.attrs["name"].decode('utf-8'), None, [], node.attrs["label"].decode('utf-8')]
            if " data" in node:
                dataitem = node[" data"]
                res[1] = np.copy(dataitem[()].transpose(), order="F")
                if node.attrs["type"] == b"C1":
                    res[1] = np.vectorize(chr)(res[1]).astype(np.dtype("c"))

            names = [x for x in node.keys() if x[0] != " "]
            for name in names:
                child = ConvertData(node[name])
                if child is not None:
                    res[2].append(child)
            return res

        node = ConvertData(h5file)
        node[0] = "CGNSTree"
        node[3] = "CGNSTree_t"
        self.CGNSTree = node

        res = CGNSToMesh(node, baseNumberOrName, zoneNumberOrName)
        res.PrepareForOutput()
        return res

from BasicTools.IO.IOFactory import RegisterReaderClass
RegisterReaderClass(".cgns", CGNSReader)

def CheckIntegrity(GUI=False):

    try:
        from h5py import File
    except:
        return "skip! h5py module not available"

    import BasicTools.IO.CGNSWriter as CW
    CW.CheckIntegrity()

    from BasicTools.Helpers.Tests import TestTempDir
    tempdir = TestTempDir.GetTempPath()

    mesh = ReadCGNS(fileName = tempdir+os.sep+"toto.cgns", baseNumberOrName = 0, zoneNumberOrName = 0)

    print("Read mesh from cgns:", mesh)

    return 'ok'

if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover
